#include "fstsp.hpp"
#include "timing.hpp"

#include "gurobi_c++.h"

#include <algorithm>
#include <cmath>
#include <exception>
#include <stdexcept>
#include <utility>

namespace fstsp {
namespace {

struct Assignment {
    int stop{}, parcel{}, variable{};
};

class BendersCallback final : public GRBCallback {
public:
    BendersCallback(const TimeGraph& graph, const std::vector<GRBVar>& variables, GRBVar eta,
                    double baseline, Stats& stats)
        : graph_(graph), variables_(variables), eta_(eta), baseline_(baseline), stats_(stats) {}

    TimeAnswer best;
    std::vector<int> best_assignment;
    std::exception_ptr error;

protected:
    void callback() override {
        if (where != GRB_CB_MIPSOL) {
            return;
        }
        try {
            std::vector<int> assignment;
            assignment.reserve(variables_.size());
            for (const auto& variable : variables_) {
                assignment.push_back(getSolution(variable) > 0.5 ? 1 : 0);
            }
            auto answer = graph_.solve(assignment, true);
            if (!answer.feasible) {
                addLazy(cut(answer.certificate, 0), GRB_LESS_EQUAL, 0);
                ++stats_.feasibility_cuts;
                return;
            }
            const double value = answer.times[graph_.finish()];
            if (!best.feasible || value < best.times[graph_.finish()]) {
                best = answer;
                best_assignment = std::move(assignment);
            }
            if (getSolution(eta_) + tolerance < value - baseline_) {
                addLazy(eta_ - cut(answer.certificate, -baseline_), GRB_GREATER_EQUAL, 0);
                ++stats_.optimality_cuts;
            }
        } catch (...) {
            error = std::current_exception();
            abort();
        }
    }

private:
    GRBLinExpr cut(const std::vector<int>& certificate, double offset) const {
        // Sum along a path for (A.39), or a positive cycle for (A.40).
        GRBLinExpr expression = offset;
        for (const int e : certificate) {
            const auto& arc = graph_.arcs[e];
            expression += arc.base;
            for (const auto& term : arc.terms) {
                if (term.complement) {
                    expression += term.coefficient * (1 - variables_[term.variable]);
                } else {
                    expression += term.coefficient * variables_[term.variable];
                }
            }
        }
        return expression;
    }

    const TimeGraph& graph_;
    const std::vector<GRBVar>& variables_;
    GRBVar eta_;
    double baseline_;
    Stats& stats_;
};

double horizon(const Instance& instance, const Structure& plan, const DeliveryInfo& info,
               const std::vector<int>& delayed) {
    double bound = 0;
    for (const int j : delayed) {
        bound = std::max(bound, instance.customers[j].release);
    }
    for (std::size_t k = 1; k < plan.truck.size(); ++k) {
        const int node = plan.truck[k];
        bound = std::max({bound, instance.road[plan.truck[k - 1]][node], instance.customers[node].truck_service});
    }
    for (const double flight : info.flight) {
        bound = std::max(bound, flight);
    }
    for (const int k : info.rendezvous) {
        const int node = plan.truck[k];
        const double tau = instance.air[0][node];
        bound = std::max({bound, 2 * tau, tau + instance.customers[node].truck_service});
    }
    return (2 * plan.truck.size() + 2 * plan.sorties.size() + info.rendezvous.size()) * bound;
}

}  // namespace

FixedResult solve_fixed(GRBEnv& env, const Instance& instance, const Structure& plan,
                       const std::vector<Pair>& warm, Deadline deadline, Stats& stats) {
    const auto info = describe(instance, plan, true);
    if (!info) {
        return {};
    }
    if (!prefix_feasible(instance, plan, *info)) {
        ++stats.screen_rejections;
        return {};
    }
    const auto baseline = delivery_schedule(instance, plan, *info, false);
    if (!baseline) {
        return {};
    }
    std::vector<int> delayed;
    for (int j = 1; j <= instance.size(); ++j) {
        if (instance.customers[j].delayed) {
            delayed.push_back(j);
        }
    }
    ++stats.fixed_solves;
    if (delayed.empty()) {
        return {FixedStatus::optimal, Solution{plan, *baseline, baseline->makespan()}};
    }
    if (Clock::now() >= deadline) {
        return {FixedStatus::time_limit, std::nullopt};
    }
    const int m = static_cast<int>(info->rendezvous.size());
    const int truck_count = static_cast<int>(plan.truck.size());
    const double time_bound = horizon(instance, plan, *info, delayed);
    const double base = baseline->makespan();
    GRBModel master(env);
    master.set(GRB_IntParam_LazyConstraints, 1);
    master.set(GRB_IntParam_Threads, 1);
    master.set(GRB_IntParam_DualReductions, 0);
    master.set(GRB_DoubleParam_MIPGap, 0);
    master.set(GRB_DoubleParam_MIPGapAbs, 0);
    master.set(GRB_DoubleParam_FeasibilityTol, 1e-8);
    master.set(GRB_DoubleParam_IntFeasTol, 1e-9);
    const auto eta = master.addVar(0, GRB_INFINITY, 1, GRB_CONTINUOUS, "eta");
    std::vector<GRBVar> variables;
    std::vector<Assignment> assignments;
    std::vector<std::vector<int>> pair_index(m, std::vector<int>(instance.customers.size(), -1));
    for (int i = 0; i < m; ++i) {
        variables.push_back(master.addVar(0, 1, 0, GRB_BINARY));
    }
    for (int i = 0; i < m; ++i) {
        for (const int j : delayed) {
            if (info->rendezvous[i] <= info->cutoff[j]) {
                const int index = static_cast<int>(variables.size());
                variables.push_back(master.addVar(0, 1, 0, GRB_BINARY));
                pair_index[i][j] = index;
                assignments.push_back({i, j, index});
            }
        }
    }
    std::vector<GRBLinExpr> count(m), weight(m), covered(instance.customers.size());
    std::vector<double> remaining(truck_count);
    for (int k = truck_count - 2; k >= 0; --k) {
        const int next = plan.truck[k + 1];
        remaining[k] = remaining[k + 1] + instance.road[plan.truck[k]][next] +
                       instance.customers[next].truck_service;
    }
    TimeGraph graph(truck_count, static_cast<int>(plan.sorties.size()), m);
    graph.delivery_rows(instance, plan, *info, false);
    for (int v = 1; v < graph.vertices; ++v) {
        graph.add(v, 0, -time_bound);
    }
    for (const auto& pair : assignments) {
        const int i = pair.stop, j = pair.parcel, v = pair.variable;
        const int position = info->rendezvous[i], node = plan.truck[position];
        const double tau = instance.air[0][node];
        const auto& parcel = instance.customers[j];
        count[i] += variables[v];
        weight[i] += parcel.weight * variables[v];
        covered[j] += variables[v];
        // (35), (60a), and (63).
        graph.add(0, graph.resupply(i), 0, {{v, parcel.release}});
        const int s = info->sortie_of[j];
        if (s >= 0) {
            graph.add(graph.resupply(i), graph.launch(s), tau, {{v, -(time_bound + tau), true}});
        }
        const double bound = parcel.release + tau + (node == j ? parcel.truck_service : 0) +
                             remaining[position] - base;
        if (bound > 0) {
            master.addConstr(eta >= bound * variables[v]);
        }
    }
    double total = 0;
    for (const int j : delayed) {
        master.addConstr(covered[j] == 1);
        total += instance.customers[j].weight;
    }
    GRBLinExpr active;
    for (int i = 0; i < m; ++i) {
        const int k = info->rendezvous[i], node = plan.truck[k];
        const double tau = instance.air[0][node];
        const int same = pair_index[i][node];
        master.addConstr(count[i] >= variables[i]);
        master.addConstr(count[i] <= static_cast<double>(delayed.size()) * variables[i]);
        master.addConstr(weight[i] <= instance.resupply_capacity * variables[i]);
        active += variables[i];

        // (50)-(52). An inactive candidate only propagates the previous return time.
        if (i + 1 < m) {
            graph.add(graph.resupply(i), graph.resupply(i + 1), 0, {{i, 2 * tau}});
        }
        std::vector<Term> dwell{{i, -(time_bound + tau), true}};
        if (same >= 0) {
            dwell.push_back({same, instance.customers[node].truck_service});
        }
        graph.add(graph.resupply(i), graph.departure(k), tau, std::move(dwell));
        graph.add(graph.arrival(k), graph.resupply(i), -tau,
                  {{i, -std::max(0.0, time_bound - tau), true}});
    }
    const double min_trips = std::max(1.0, std::ceil(total / instance.resupply_capacity - 1e-10));
    master.addConstr(active >= min_trips);
    master.addConstr(active <= static_cast<double>(std::min(info->rendezvous.size(), delayed.size())));
    master.update();
    int start_count = 0;
    for (const auto& pair : warm) {
        const int k = info->position[pair.rendezvous];
        const auto found = std::lower_bound(info->rendezvous.begin(), info->rendezvous.end(), k);
        if (found != info->rendezvous.end() && *found == k) {
            const int i = static_cast<int>(found - info->rendezvous.begin());
            const int index = pair_index[i][pair.parcel];
            if (index >= 0) {
                variables[index].set(GRB_DoubleAttr_Start, 1);
                variables[i].set(GRB_DoubleAttr_Start, 1);
                ++start_count;
            }
        }
    }
    stats.warm_starts += start_count > 0;
    stats.warm_pairs += start_count;
    const double remaining_seconds = std::chrono::duration<double>(deadline - Clock::now()).count();
    if (remaining_seconds <= 0) {
        return {FixedStatus::time_limit, std::nullopt};
    }
    master.set(GRB_DoubleParam_TimeLimit, remaining_seconds);
    BendersCallback callback(graph, variables, eta, base, stats);
    master.setCallback(&callback);
    master.optimize();
    master.setCallback(nullptr);
    if (callback.error) {
        std::rethrow_exception(callback.error);
    }
    const int status = master.get(GRB_IntAttr_Status);
    if (status == GRB_INFEASIBLE) {
        if (callback.best.feasible) {
            throw std::runtime_error("Benders master excluded a feasible timetable.");
        }
        return {};
    }
    if (status != GRB_OPTIMAL && status != GRB_TIME_LIMIT) {
        throw std::runtime_error("Unexpected fixed-structure solver status: " + std::to_string(status));
    }
    if (!callback.best.feasible) {
        if (status == GRB_OPTIMAL) {
            throw std::runtime_error("An optimal master solution was not scheduled.");
        }
        return {FixedStatus::time_limit, std::nullopt};
    }
    const auto& times = callback.best.times;
    if (status == GRB_OPTIMAL &&
        std::abs(master.get(GRB_DoubleAttr_ObjVal) - (times[graph.finish()] - base)) > 1e-5) {
        throw std::runtime_error("Benders master and exact timing objective disagree.");
    }
    auto schedule = graph.timetable(times);
    std::vector<std::vector<int>> parcels(m);
    for (const auto& assignment : assignments) {
        if (callback.best_assignment[assignment.variable]) {
            parcels[assignment.stop].push_back(assignment.parcel);
        }
    }
    for (int i = 0; i < m; ++i) {
        if (!parcels[i].empty()) {
            const int node = plan.truck[info->rendezvous[i]];
            const double depart = times[graph.resupply(i)], tau = instance.air[0][node];
            schedule.resupply.push_back({node, depart, depart + tau, depart + 2 * tau, std::move(parcels[i])});
        }
    }
    Solution answer{plan, std::move(schedule), base};
    verify_solution(instance, answer);
    return {status == GRB_OPTIMAL ? FixedStatus::optimal : FixedStatus::time_limit, std::move(answer)};
}

}  // namespace fstsp
