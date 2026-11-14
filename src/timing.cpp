#include "timing.hpp"

#include <algorithm>
#include <stdexcept>
#include <utility>

namespace fstsp {

TimeGraph::TimeGraph(int size, int sorties, int stops)
    : truck_count(size), sortie_count(sorties), vertices(1 + 2 * size + 2 * sorties + stops) {
    for (int i = 1; i < vertices; ++i) {
        add(0, i, 0);
    }
}

void TimeGraph::add(int from, int to, double rhs, std::vector<Term> terms) {
    arcs.push_back({from, to, rhs, std::move(terms)});
}

void TimeGraph::delivery_rows(const Instance& instance, const Structure& plan, const DeliveryInfo& info,
                              bool availability_bounds) {
    add(arrival(0), 0, 0);
    for (int k = 1; k < truck_count; ++k) {
        const int previous = plan.truck[k - 1], node = plan.truck[k];
        const double drive = instance.road[previous][node];
        const auto& customer = instance.customers[node];
        add(departure(k - 1), arrival(k), drive);
        add(arrival(k), departure(k - 1), -drive);
        add(arrival(k), departure(k), customer.truck_service);
        if (availability_bounds && customer.delayed) {
            add(0, departure(k), customer.release + instance.air[0][node] + customer.truck_service);
        }
    }
    for (int s = 0; s < sortie_count; ++s) {
        const int l = info.launch_pos[s], b = info.recovery_pos[s];
        add(arrival(l), launch(s), 0);
        add(launch(s), departure(l), 0);
        add(arrival(b), recovery(s), 0);
        add(recovery(s), departure(b), 0);
        add(launch(s), recovery(s), info.flight[s]);
        add(recovery(s), launch(s), -instance.delivery_endurance);
        if (s > 0) {
            add(recovery(s - 1), launch(s), 0);
        }
        if (availability_bounds) {
            for (const int j : plan.sorties[s].visits) {
                if (instance.customers[j].delayed) {
                    add(0, launch(s), instance.customers[j].release + instance.air[0][plan.truck[l]]);
                }
            }
        }
    }
}

TimeAnswer TimeGraph::solve(const std::vector<int>& binary, bool get_certificate) const {
    constexpr double epsilon = 1e-9;
    std::vector<double> value(vertices), rhs;
    std::vector<int> previous(vertices, -1);
    rhs.reserve(arcs.size());
    for (const auto& arc : arcs) {
        double weight = arc.base;
        for (const auto& term : arc.terms) {
            if ((binary[term.variable] != 0) != term.complement) {
                weight += term.coefficient;
            }
        }
        rhs.push_back(weight);
    }
    int changed = -1;
    for (int pass = 0; pass < vertices; ++pass) {
        changed = -1;
        for (int e = 0; e < static_cast<int>(arcs.size()); ++e) {
            const auto& arc = arcs[e];
            const double candidate = value[arc.from] + rhs[e];
            if (candidate > value[arc.to] + epsilon) {
                value[arc.to] = candidate;
                previous[arc.to] = e;
                changed = arc.to;
            }
        }
        if (changed < 0) {
            break;
        }
    }
    if (!get_certificate) {
        return changed < 0 ? TimeAnswer{true, std::move(value), {}} : TimeAnswer{};
    }
    if (changed < 0) {
        // The predecessor chain gives an optimal dual path from the fixed time zero.
        std::vector<int> path;
        for (int v = finish(); v != 0;) {
            const int edge = previous[v];
            if (edge < 0) {
                path.push_back(v - 1);  // The nonnegativity arc 0 -> v.
                break;
            }
            if (path.size() >= static_cast<std::size_t>(vertices)) {
                throw std::runtime_error("Could not recover a timing path.");
            }
            path.push_back(edge);
            v = arcs[edge].from;
        }
        return {true, std::move(value), std::move(path)};
    }

    // An improving relaxation on pass |V| implies a positive cycle.
    int cycle_vertex = changed;
    for (int step = 0; step < vertices; ++step) {
        if (previous[cycle_vertex] < 0) {
            throw std::runtime_error("Could not recover a timing infeasibility certificate.");
        }
        cycle_vertex = arcs[previous[cycle_vertex]].from;
    }
    std::vector<int> cycle;
    int vertex = cycle_vertex;
    do {
        const int edge = previous[vertex];
        cycle.push_back(edge);
        vertex = arcs[edge].from;
    } while (vertex != cycle_vertex);
    return {false, {}, std::move(cycle)};
}

Schedule TimeGraph::timetable(const std::vector<double>& times) const {
    Schedule result;
    for (int k = 0; k < truck_count; ++k) {
        result.arrival.push_back(times[arrival(k)]);
        result.departure.push_back(times[departure(k)]);
    }
    for (int s = 0; s < sortie_count; ++s) {
        result.launch.push_back(times[launch(s)]);
        result.recovery.push_back(times[recovery(s)]);
    }
    return result;
}

std::optional<Schedule> delivery_schedule(const Instance& instance, const Structure& plan,
                                          const DeliveryInfo& info, bool availability_bounds) {
    TimeGraph graph(static_cast<int>(plan.truck.size()), static_cast<int>(plan.sorties.size()));
    graph.delivery_rows(instance, plan, info, availability_bounds);
    const auto result = graph.solve({});
    if (!result.feasible) {
        return std::nullopt;
    }
    return graph.timetable(result.times);
}

}  // namespace fstsp
