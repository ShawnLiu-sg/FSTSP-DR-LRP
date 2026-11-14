#include "fstsp.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <random>
#include <utility>

namespace fstsp {
namespace {

using Random = std::mt19937_64;
enum class Insertion { truck, new_sortie, extend_sortie };
enum class InsertionRule { greedy, noisy, random_customer };
enum class Removal { random, shaw, worst };

struct Move {
    Insertion type;
    int node, first, second;
    double raw, score;
};

struct Tentative {
    Structure plan;
    DeliveryInfo info;
    Schedule times;
    std::vector<double> wait, suffix;
    std::vector<int> covered_edge;

    void update_suffix() {
        suffix.assign(wait.size(), 0);
        for (int k = static_cast<int>(wait.size()) - 2; k >= 0; --k) {
            suffix[k] = suffix[k + 1] + wait[k + 1];
        }
    }
};

struct Placed {
    Tentative state;
    std::vector<int> displaced;
    bool changed_recovery{};
};

class Operators {
public:
    Operators(int count, int segment, double reaction)
        : weights_(count, 1), scores_(count), uses_(count), segment_(segment), reaction_(reaction) {}

    int choose(Random& random) const {
        std::discrete_distribution<int> distribution(weights_.begin(), weights_.end());
        return distribution(random);
    }

    void reward(int which, double score) {
        scores_[which] += score;
        ++uses_[which];
        if (++steps_ == segment_) {
            for (std::size_t i = 0; i < weights_.size(); ++i) {
                if (uses_[i] > 0) {
                    weights_[i] += reaction_ * (scores_[i] / uses_[i] - weights_[i]);
                }
            }
            std::fill(scores_.begin(), scores_.end(), 0);
            std::fill(uses_.begin(), uses_.end(), 0);
            steps_ = 0;
        }
    }

private:
    std::vector<double> weights_, scores_;
    std::vector<int> uses_;
    int segment_;
    double reaction_;
    int steps_{};
};

double random_unit(Random& random) {
    return std::generate_canonical<double, 53>(random);
}

bool retained(const DeliveryInfo& info, const std::vector<Pair>& keep) {
    for (const auto& pair : keep) {
        const int rendezvous = info.position[pair.rendezvous], cutoff = info.cutoff[pair.parcel];
        if (info.sortie_of[pair.rendezvous] >= 0 ||
            (rendezvous >= 0 && cutoff >= 0 && rendezvous > cutoff)) {
            return false;
        }
    }
    return true;
}

std::optional<Tentative> tentative(const Instance& instance, const Structure& plan,
                                   const std::vector<Pair>& keep, bool force_screen) {
    auto info = describe(instance, plan, false);
    if (!info || !retained(*info, keep) || (force_screen && !prefix_feasible(instance, plan, *info))) {
        return std::nullopt;
    }
    auto times = delivery_schedule(instance, plan, *info, true);
    if (!times) {
        return std::nullopt;
    }
    Tentative result{plan, std::move(*info), std::move(*times), {}, {}, {}};
    const auto& t = result.times;
    std::vector<double> ready = t.arrival;
    for (std::size_t k = 0; k < ready.size(); ++k) {
        ready[k] += instance.customers[plan.truck[k]].truck_service;
    }
    result.covered_edge.assign(plan.truck.size() - 1, -1);
    for (int s = 0; s < static_cast<int>(plan.sorties.size()); ++s) {
        const int l = result.info.launch_pos[s], b = result.info.recovery_pos[s];
        ready[l] = std::max(ready[l], t.recovery[s] - instance.delivery_endurance);
        ready[b] = std::max(ready[b], t.recovery[s]);
        std::fill(result.covered_edge.begin() + l, result.covered_edge.begin() + b, s);
    }
    for (std::size_t k = 0; k < ready.size(); ++k) {
        result.wait.push_back(std::max(0.0, t.departure[k] - ready[k]));
    }
    result.update_suffix();
    return result;
}

void propagate_slack(const Tentative& old, const Move& move, Placed& placed) {
    if (placed.changed_recovery || !placed.displaced.empty()) {
        return;  // A new recovery structure has its own tentative stop windows.
    }
    const int begin = move.type == Insertion::extend_sortie ? old.info.launch_pos[move.first] : move.first;
    double residual = std::max(0.0, move.raw);
    for (int k = begin + 1; k < static_cast<int>(old.wait.size()); ++k) {
        const int next = placed.state.info.position[old.plan.truck[k]];
        const double absorbed = std::min(residual, old.wait[k]);
        placed.state.wait[next] = old.wait[k] - absorbed;
        residual -= absorbed;
    }
    placed.state.update_suffix();
}

double launch_bound(const Instance& instance, int node, int launch) {
    return instance.customers[node].delayed ? instance.customers[node].release + instance.air[0][launch] : 0;
}

void insertion_moves(const Instance& instance, const Tentative& state, int node, bool truck_only,
                     double noise, Random& random, std::vector<Move>& moves) {
    const auto& customer = instance.customers[node];
    const auto& plan = state.plan;
    const auto& t = state.times;
    const auto& info = state.info;
    const int count = static_cast<int>(plan.truck.size());
    const auto add = [&](Insertion type, int first, int second, int k, double raw) {
        const double net = std::max(0.0, raw - state.suffix[k]);
        const double score = net + (noise > 0 ? net * noise * (2 * random_unit(random) - 1) : 0);
        moves.push_back({type, node, first, second, raw, score});
    };
    for (int j = 0; j + 1 < count; ++j) {
        const int a = plan.truck[j], b = plan.truck[j + 1];
        const double raw = instance.road[a][node] + customer.truck_service +
                           instance.road[node][b] - instance.road[a][b] +
                           std::max(0.0, launch_bound(instance, node, node) -
                                         t.departure[j] - instance.road[a][node]);
        add(Insertion::truck, j, state.covered_edge[j], j, raw);
    }
    if (truck_only || customer.weight > instance.delivery_capacity + tolerance) {
        return;
    }
    for (int j = 0; j + 1 < count; ++j) {
        const int a = plan.truck[j];
        const double before = std::max(0.0, launch_bound(instance, node, a) - t.departure[j]);
        for (int k = j + 1; k < count && state.covered_edge[k - 1] < 0; ++k) {
            const int b = plan.truck[k];
            const double duration = instance.air[a][node] + customer.drone_service + instance.air[node][b];
            if (duration <= instance.delivery_endurance + tolerance) {
                const double window = t.arrival[k] - t.departure[j];
                add(Insertion::new_sortie, j, k, j, before + std::max(0.0, duration - window));
            }
        }
    }
    for (int s = 0; s < static_cast<int>(plan.sorties.size()); ++s) {
        if (info.payload[s] + customer.weight > instance.delivery_capacity + tolerance) {
            continue;
        }
        const auto& sortie = plan.sorties[s];
        const double before = std::max(0.0, launch_bound(instance, node, sortie.launch) - t.launch[s]);
        const double window = t.arrival[info.recovery_pos[s]] - t.launch[s];
        int a = sortie.launch;
        for (int pos = 0; pos <= static_cast<int>(sortie.visits.size()); ++pos) {
            const int b = pos < static_cast<int>(sortie.visits.size()) ? sortie.visits[pos] : sortie.recovery;
            const double duration = info.flight[s] + instance.air[a][node] + customer.drone_service +
                                    instance.air[node][b] - instance.air[a][b];
            if (duration <= instance.delivery_endurance + tolerance) {
                const double raw = before + std::max(0.0, duration - window) -
                                   std::max(0.0, info.flight[s] - window);
                add(Insertion::extend_sortie, s, pos, info.launch_pos[s], raw);
            }
            a = b;
        }
    }
}

std::pair<int, int> recovery_range(const Structure& plan, int s) {
    const auto& sortie = plan.sorties[s];
    const auto launch = std::find(plan.truck.begin(), plan.truck.end(), sortie.launch);
    const auto latest = s + 1 < static_cast<int>(plan.sorties.size())
        ? std::find(plan.truck.begin(), plan.truck.end(), plan.sorties[s + 1].launch)
        : plan.truck.end() - 1;
    return {static_cast<int>(launch - plan.truck.begin()) + 1, static_cast<int>(latest - plan.truck.begin())};
}

std::optional<Placed> place(const Instance& instance, const Tentative& old, const Move& move,
                            const std::vector<Pair>& keep, bool force_screen, Deadline deadline) {
    Structure trial = old.plan;
    if (move.type == Insertion::new_sortie) {
        const int launch = trial.truck[move.first], recovery = trial.truck[move.second];
        const auto first = std::lower_bound(old.info.launch_pos.begin(), old.info.launch_pos.end(), move.first);
        trial.sorties.insert(trial.sorties.begin() + (first - old.info.launch_pos.begin()),
                             {launch, recovery, {move.node}});
    } else if (move.type == Insertion::extend_sortie) {
        auto& visits = trial.sorties[move.first].visits;
        visits.insert(visits.begin() + move.second, move.node);
    } else {
        trial.truck.insert(trial.truck.begin() + move.first + 1, move.node);
    }
    if (move.type != Insertion::truck || move.second < 0) {
        if (auto state = tentative(instance, trial, keep, force_screen)) {
            return Placed{std::move(*state), {}, false};
        }
        return std::nullopt;
    }

    const int s = move.second;
    const int original = trial.sorties[s].recovery;
    const auto [first, last] = recovery_range(trial, s);
    std::optional<Placed> best;
    double best_duration = std::numeric_limits<double>::infinity();
    for (int k = first; k <= last && Clock::now() < deadline; ++k) {
        trial.sorties[s].recovery = trial.truck[k];
        auto state = tentative(instance, trial, keep, force_screen);
        if (state) {
            const auto baseline = delivery_schedule(instance, trial, state->info, false);
            if (baseline->makespan() < best_duration - tolerance ||
                (trial.truck[k] == original && baseline->makespan() <= best_duration + tolerance)) {
                best_duration = baseline->makespan();
                best = Placed{std::move(*state), {}, trial.truck[k] != original};
            }
        }
    }
    if (best) {
        return best;
    }
    auto displaced = std::move(trial.sorties[s].visits);
    trial.sorties.erase(trial.sorties.begin() + s);
    if (auto state = tentative(instance, trial, keep, false)) {
        return Placed{std::move(*state), std::move(displaced), true};
    }
    return std::nullopt;
}

std::optional<Structure> insert(const Instance& instance, Structure partial, std::vector<int> pending,
                                InsertionRule rule, const std::vector<Pair>& keep, bool force_screen,
                                const Settings& settings, Random& random, Deadline deadline) {
    auto current = tentative(instance, partial, keep, false);
    if (!current) {
        return std::nullopt;
    }
    std::vector<bool> truck_only(instance.customers.size());
    for (const auto& pair : keep) {
        truck_only[pair.rendezvous] = true;
    }
    const auto try_moves = [&](std::vector<Move>& moves) -> std::optional<Placed> {
        std::stable_sort(moves.begin(), moves.end(), [](const Move& a, const Move& b) { return a.score < b.score; });
        for (const auto& move : moves) {
            if (Clock::now() >= deadline) {
                break;
            }
            auto result = place(instance, *current, move, keep, force_screen && pending.size() == 1, deadline);
            if (result) {
                propagate_slack(*current, move, *result);
                pending.erase(std::find(pending.begin(), pending.end(), move.node));
                return result;
            }
        }
        return std::nullopt;
    };
    while (!pending.empty() && Clock::now() < deadline) {
        const bool need_truck = current->plan.truck.size() == 2 && pending.size() == 1;
        std::optional<Placed> next;
        if (rule == InsertionRule::random_customer) {
            auto order = pending;
            std::shuffle(order.begin(), order.end(), random);
            for (const int node : order) {
                std::vector<Move> moves;
                insertion_moves(instance, *current, node, need_truck || truck_only[node], 0, random, moves);
                next = try_moves(moves);
                if (next || Clock::now() >= deadline) {
                    break;
                }
            }
        } else {
            std::vector<Move> moves;
            for (const int node : pending) {
                if (Clock::now() >= deadline) {
                    return std::nullopt;
                }
                const double noise = rule == InsertionRule::noisy ? settings.noise : 0;
                insertion_moves(instance, *current, node, need_truck || truck_only[node], noise, random, moves);
            }
            next = try_moves(moves);
        }
        if (!next) {
            return std::nullopt;
        }
        for (const int node : next->displaced) {
            truck_only[node] = true;
            pending.push_back(node);
        }
        current = std::move(next->state);
    }
    if (!pending.empty() || current->plan.truck.size() < 3 ||
        (force_screen && !prefix_feasible(instance, current->plan, current->info))) {
        return std::nullopt;
    }
    return std::move(current->plan);
}

std::optional<Schedule> route_timing(const Instance& instance, const Structure& plan, bool complete) {
    const auto info = describe(instance, plan, complete);
    return info ? delivery_schedule(instance, plan, *info, false) : std::nullopt;
}

double delivery_value(const Instance& instance, const Structure& plan) {
    const auto schedule = route_timing(instance, plan, false);
    return schedule ? schedule->makespan() : std::numeric_limits<double>::infinity();
}

std::vector<int> customers(const Instance& instance) {
    std::vector<int> result(instance.size());
    std::iota(result.begin(), result.end(), 1);
    return result;
}

std::vector<std::vector<int>> cluster_customers(const Instance& instance, const Settings& settings,
                                                Random& random, Deadline deadline) {
    const auto nodes = customers(instance);
    std::vector<std::vector<int>> best{nodes};
    double best_silhouette = 0;
    const int largest_k = std::min(settings.max_clusters, instance.size() - 1);
    using Point = std::pair<double, double>;
    const auto distance_squared = [](const Point& a, const Point& b) {
        const double dx = a.first - b.first, dy = a.second - b.second;
        return dx * dx + dy * dy;
    };
    std::vector<Point> points;
    for (const int node : nodes) {
        const auto& customer = instance.customers[node];
        points.emplace_back(customer.x, customer.y);
    }
    for (int k = 2; k <= largest_k && Clock::now() < deadline; ++k) {
        std::uniform_int_distribution<int> first(0, instance.size() - 1);
        std::vector<Point> centers{points[first(random)]};
        std::vector<double> distances(points.size(), std::numeric_limits<double>::infinity());
        while (static_cast<int>(centers.size()) < k) {
            for (std::size_t i = 0; i < points.size(); ++i) {
                distances[i] = std::min(distances[i], distance_squared(points[i], centers.back()));
            }
            if (std::all_of(distances.begin(), distances.end(), [](double d) { return d == 0; })) {
                break;
            }
            std::discrete_distribution<std::size_t> choice(distances.begin(), distances.end());
            centers.push_back(points[choice(random)]);
        }
        if (static_cast<int>(centers.size()) != k) {
            break;
        }
        std::vector<std::vector<int>> groups(k);
        std::vector<int> membership(points.size(), -1);
        for (int pass = 0; pass < settings.cluster_iterations; ++pass) {
            for (auto& group : groups) {
                group.clear();
            }
            bool changed = false;
            for (std::size_t i = 0; i < points.size(); ++i) {
                int closest = 0;
                for (int c = 1; c < k; ++c) {
                    if (distance_squared(points[i], centers[c]) < distance_squared(points[i], centers[closest])) {
                        closest = c;
                    }
                }
                changed |= membership[i] != closest;
                membership[i] = closest;
                groups[closest].push_back(nodes[i]);
            }
            if (!changed || Clock::now() >= deadline) {
                break;
            }
            for (int c = 0; c < k; ++c) {
                if (!groups[c].empty()) {
                    Point center{};
                    for (const int node : groups[c]) {
                        center.first += instance.customers[node].x;
                        center.second += instance.customers[node].y;
                    }
                    const double count = static_cast<double>(groups[c].size());
                    centers[c] = {center.first / count, center.second / count};
                }
            }
        }
        if (std::any_of(groups.begin(), groups.end(), [](const auto& group) { return group.empty(); })) {
            continue;
        }
        double silhouette = 0;
        for (std::size_t i = 0; i < nodes.size(); ++i) {
            const int own = membership[i];
            if (groups[own].size() <= 1) {
                continue;
            }
            double inside = 0, outside = std::numeric_limits<double>::infinity();
            for (int c = 0; c < k; ++c) {
                double sum = 0;
                for (const int node : groups[c]) {
                    sum += instance.air[nodes[i]][node];
                }
                if (c == own) {
                    inside = sum / (groups[c].size() - 1);
                } else {
                    outside = std::min(outside, sum / groups[c].size());
                }
            }
            if (std::max(inside, outside) > 0) {
                silhouette += (outside - inside) / std::max(inside, outside);
            }
        }
        silhouette /= static_cast<double>(nodes.size());
        if (silhouette > best_silhouette) {
            best_silhouette = silhouette;
            best = std::move(groups);
        }
    }
    return best;
}

std::optional<Structure> construct(const Instance& instance, int which, const Settings& settings,
                                   Random& random, Deadline deadline) {
    const Structure empty{{0, instance.end()}, {}};
    if (which < 2) {
        return insert(instance, empty, customers(instance),
                      which == 0 ? InsertionRule::noisy : InsertionRule::random_customer,
                      {}, false, settings, random, deadline);
    }

    auto groups = cluster_customers(instance, settings, random, deadline);
    std::vector<Structure> routes;
    for (const auto& group : groups) {
        auto route = insert(instance, empty, group, InsertionRule::noisy, {}, false, settings, random, deadline);
        if (!route) {
            return std::nullopt;
        }
        routes.push_back(std::move(*route));
    }
    Structure merged = empty;
    while (!routes.empty() && Clock::now() < deadline) {
        const auto current = tentative(instance, merged, {}, false);
        std::optional<Structure> best;
        double best_value = std::numeric_limits<double>::infinity();
        std::size_t selected = 0;
        const int edges = static_cast<int>(current->covered_edge.size());
        const int first_edge = settings.cluster_join == ClusterJoin::append ? edges - 1 : 0;
        for (std::size_t c = 0; c < routes.size() && Clock::now() < deadline; ++c) {
            const auto& route = routes[c];
            for (int k = first_edge; k < edges; ++k) {
                if (current->covered_edge[k] >= 0) {
                    continue;
                }
                Structure candidate = merged;
                candidate.truck.insert(candidate.truck.begin() + k + 1,
                                       route.truck.begin() + 1, route.truck.end() - 1);
                auto sorties = route.sorties;
                for (auto& sortie : sorties) {
                    if (sortie.launch == 0) {
                        sortie.launch = merged.truck[k];
                    }
                    if (sortie.recovery == instance.end()) {
                        sortie.recovery = merged.truck[k + 1];
                    }
                }
                const auto at = std::lower_bound(current->info.launch_pos.begin(),
                                                 current->info.launch_pos.end(), k);
                candidate.sorties.insert(candidate.sorties.begin() + (at - current->info.launch_pos.begin()),
                                         sorties.begin(), sorties.end());
                auto state = tentative(instance, candidate, {}, false);
                if (state && state->times.makespan() < best_value - tolerance) {
                    best_value = state->times.makespan();
                    best = std::move(candidate);
                    selected = c;
                }
            }
        }
        if (!best) {
            // Reinsert the cluster's customers when its subroute has no feasible splice.
            auto partial = describe(instance, merged, false);
            std::vector<int> pending;
            for (const int j : customers(instance)) {
                if (partial->cutoff[j] < 0) {
                    pending.push_back(j);
                }
            }
            return insert(instance, merged, std::move(pending), InsertionRule::greedy,
                          {}, false, settings, random, deadline);
        }
        merged = std::move(*best);
        routes.erase(routes.begin() + selected);
    }
    return routes.empty() ? std::optional<Structure>(std::move(merged)) : std::nullopt;
}

Structure remove_customers(const Structure& plan, const std::vector<bool>& removed) {
    Structure partial;
    for (const int node : plan.truck) {
        if (!removed[node]) {
            partial.truck.push_back(node);
        }
    }
    for (const auto& sortie : plan.sorties) {
        if (removed[sortie.launch] || removed[sortie.recovery]) {
            continue;  // These customers rejoin the insertion queue.
        }
        Sortie kept{sortie.launch, sortie.recovery, {}};
        for (const int node : sortie.visits) {
            if (!removed[node]) {
                kept.visits.push_back(node);
            }
        }
        if (!kept.visits.empty()) {
            partial.sorties.push_back(std::move(kept));
        }
    }
    return partial;
}

Structure destroy(const Instance& instance, const Structure& plan, Removal rule, const Settings& settings,
                  Random& random, Deadline deadline) {
    auto available = customers(instance);
    const int count = static_cast<int>(std::ceil(settings.remove_fraction * instance.size()));
    std::vector<bool> removed(instance.customers.size());
    if (rule == Removal::random) {
        std::shuffle(available.begin(), available.end(), random);
        for (int k = 0; k < count; ++k) {
            removed[available[k]] = true;
        }
    } else if (rule == Removal::shaw) {
        double distance_scale = 0, release_scale = 0, weight_scale = 0;
        for (const int i : available) {
            distance_scale = std::max(distance_scale, *std::max_element(instance.air[i].begin(), instance.air[i].end()));
            release_scale = std::max(release_scale, instance.customers[i].release);
            weight_scale = std::max(weight_scale, instance.customers[i].weight);
        }
        const auto related = [&](int a, int b) {
            const auto& p = instance.customers[a];
            const auto& q = instance.customers[b];
            return (distance_scale > 0 ? instance.air[a][b] / distance_scale : 0) +
                   (release_scale > 0 ? std::abs(p.release - q.release) / release_scale : 0) +
                   std::abs(p.weight - q.weight) / weight_scale;
        };
        std::vector<int> selected;
        for (int k = 0; k < count; ++k) {
            int index;
            if (selected.empty()) {
                std::uniform_int_distribution<int> pick(0, static_cast<int>(available.size()) - 1);
                index = pick(random);
            } else {
                std::uniform_int_distribution<int> pick(0, static_cast<int>(selected.size()) - 1);
                const int reference = selected[pick(random)];
                std::stable_sort(available.begin(), available.end(), [&](int a, int b) {
                    return related(a, reference) < related(b, reference);
                });
                index = static_cast<int>(std::pow(random_unit(random), settings.shaw_exponent) * available.size());
            }
            selected.push_back(available[index]);
            removed[available[index]] = true;
            available.erase(available.begin() + index);
        }
    } else {
        for (int k = 0; k < count && !available.empty() && Clock::now() < deadline; ++k) {
            int worst = 0;
            double best_value = std::numeric_limits<double>::infinity();
            for (int i = 0; i < static_cast<int>(available.size()) && Clock::now() < deadline; ++i) {
                removed[available[i]] = true;
                const double value = delivery_value(instance, remove_customers(plan, removed));
                removed[available[i]] = false;
                if (value < best_value - tolerance) {
                    best_value = value;
                    worst = i;
                }
            }
            removed[available[worst]] = true;
            available.erase(available.begin() + worst);
        }
    }
    return remove_customers(plan, removed);
}

std::optional<Structure> repair(const Instance& instance, Structure partial, InsertionRule rule,
                                const std::vector<Pair>& keep, bool force_screen, const Settings& settings,
                                Random& random, Deadline deadline) {
    const auto info = describe(instance, partial, false);
    if (!info) {
        return std::nullopt;
    }
    std::vector<int> pending;
    for (int j = 1; j <= instance.size(); ++j) {
        if (info->cutoff[j] < 0) {
            pending.push_back(j);
        }
    }
    return insert(instance, std::move(partial), std::move(pending), rule, keep, force_screen,
                  settings, random, deadline);
}

std::optional<Structure> neighbor(const Instance& instance, const Structure& current, int which,
                                  const std::vector<Pair>& keep, bool force_screen, const Settings& settings,
                                  Random& random, Deadline deadline, Stats& stats) {
    ++stats.iterations;
    ++stats.operator_uses[which];
    const auto removal = static_cast<Removal>(which / 2);
    const auto insertion = which % 2 == 0 ? InsertionRule::greedy : InsertionRule::random_customer;
    return repair(instance, destroy(instance, current, removal, settings, random, deadline),
                  insertion, keep, force_screen, settings, random, deadline);
}

enum class MoveOutcome { rejected, accepted, best };

MoveOutcome accept_move(double value, double current, double best, double temperature, Random& random) {
    if (!std::isfinite(value)) {
        return MoveOutcome::rejected;
    }
    if (value < best - tolerance) {
        return MoveOutcome::best;
    }
    if (value <= current + tolerance ||
        (temperature > 0 && random_unit(random) < std::exp((current - value) / temperature))) {
        return MoveOutcome::accepted;
    }
    return MoveOutcome::rejected;
}

double move_score(MoveOutcome outcome, const Settings& settings) {
    switch (outcome) {
        case MoveOutcome::best: return settings.local_best_score;
        case MoveOutcome::accepted: return settings.local_accepted_score;
        case MoveOutcome::rejected: return settings.local_rejected_score;
    }
    return settings.local_rejected_score;
}

double initial_temperature(double value, const Settings& settings) {
    return settings.initial_worsening * value / -std::log(settings.initial_acceptance);
}

std::optional<Solution> complete(GRBEnv& env, const Instance& instance, const Structure& plan,
                                const std::vector<Pair>& keep, const Settings& settings, Random& random,
                                Deadline deadline, Stats& stats) {
    auto fixed = solve_fixed(env, instance, plan, keep, deadline, stats);
    if (fixed.solution || fixed.status == FixedStatus::time_limit || plan.sorties.empty()) {
        return std::move(fixed.solution);
    }

    // A.4: each recovery choice has its own fixed delivery structure and its own cuts.
    std::optional<Solution> best;
    const auto consider = [&](const Structure& candidate) {
        const auto info = describe(instance, candidate, true);
        if (!info || !retained(*info, keep) || !prefix_feasible(instance, candidate, *info)) {
            return;
        }
        ++stats.recovery_trials;
        auto result = solve_fixed(env, instance, candidate, keep, deadline, stats);
        if (result.solution && (!best || result.solution->timing.makespan() < best->timing.makespan() - tolerance)) {
            best = std::move(result.solution);
        }
    };
    for (int s = 0; s < static_cast<int>(plan.sorties.size()) && Clock::now() < deadline; ++s) {
        const auto [first, last] = recovery_range(plan, s);
        Structure trial = plan;
        for (int k = first; k <= last && Clock::now() < deadline; ++k) {
            if (plan.truck[k] != plan.sorties[s].recovery) {
                trial.sorties[s].recovery = plan.truck[k];
                consider(trial);
            }
        }
    }
    if (best) {
        return best;
    }
    for (int s = 0; s < static_cast<int>(plan.sorties.size()) && Clock::now() < deadline; ++s) {
        Structure partial = plan;
        partial.sorties.erase(partial.sorties.begin() + s);
        if (auto trial = repair(instance, std::move(partial), InsertionRule::greedy, keep, true,
                                settings, random, deadline)) {
            consider(*trial);
        }
    }
    return best;
}

Solution local_search(GRBEnv& env, const Instance& instance, const Solution& initial,
                      const Settings& settings, Random& random, Deadline deadline, Stats& stats) {
    Solution reference = initial, best = initial;
    Structure current = initial.delivery;
    double best_delivery = initial.baseline, current_delivery = initial.baseline;
    const double p0 = 1.0 / instance.size();
    double p = p0;
    double temperature = initial_temperature(initial.baseline, settings);
    int outside = 0;
    bool force = false;
    Operators operators(6, settings.weight_segment, settings.weight_reaction);
    auto permutation = assignments(reference);
    std::shuffle(permutation.begin(), permutation.end(), random);
    while (p < settings.max_relax - tolerance && Clock::now() < deadline) {
        const int which = operators.choose(random);
        const auto count = static_cast<std::size_t>(std::ceil((1 - p) * permutation.size() - 1e-12));
        const std::vector<Pair> keep(permutation.begin(), permutation.begin() + count);
        const bool force_now = std::exchange(force, false);
        stats.forced_repairs += force_now;
        const auto candidate = neighbor(instance, current, which, keep, force_now,
                                        settings, random, deadline, stats);
        bool improved = false;
        auto outcome = MoveOutcome::rejected;
        if (candidate) {
            const double value = delivery_value(instance, *candidate);
            outcome = accept_move(value, current_delivery, best_delivery, temperature, random);
            if (outcome != MoveOutcome::rejected) {
                if (outcome == MoveOutcome::best) {
                    best_delivery = value;
                }
                const auto info = describe(instance, *candidate, true);
                if (!prefix_feasible(instance, *candidate, *info)) {
                    current = *candidate;
                    current_delivery = value;
                    ++stats.screen_moves;
                    ++stats.screen_rejections;
                    if (++outside >= settings.max_infeasible_moves) {
                        current = reference.delivery;
                        current_delivery = reference.baseline;
                        outside = 0;
                        force = true;
                    }
                } else {
                    auto result = complete(env, instance, *candidate, keep, settings, random, deadline, stats);
                    if (result) {
                        const auto next_assignments = assignments(*result);
                        const auto old_assignments = assignments(reference);
                        // Keep the same permutation while the complete reference is unchanged.
                        const auto equal_assignments = [](const Pair& a, const Pair& b) {
                            return a.rendezvous == b.rendezvous && a.parcel == b.parcel;
                        };
                        const auto equal_sorties = [](const Sortie& a, const Sortie& b) {
                            return a.launch == b.launch && a.recovery == b.recovery && a.visits == b.visits;
                        };
                        const bool new_reference = reference.delivery.truck != result->delivery.truck ||
                            !std::equal(reference.delivery.sorties.begin(), reference.delivery.sorties.end(),
                                        result->delivery.sorties.begin(), result->delivery.sorties.end(), equal_sorties) ||
                            !std::equal(old_assignments.begin(), old_assignments.end(),
                                        next_assignments.begin(), next_assignments.end(), equal_assignments);
                        reference = std::move(*result);
                        current = reference.delivery;
                        current_delivery = reference.baseline;
                        best_delivery = std::min(best_delivery, current_delivery);
                        outside = 0;
                        if (new_reference) {
                            permutation = next_assignments;
                            std::shuffle(permutation.begin(), permutation.end(), random);
                        }
                        if (reference.timing.makespan() < best.timing.makespan() - tolerance) {
                            best = reference;
                            improved = true;
                            ++stats.complete_improvements;
                        }
                    } else {
                        current = reference.delivery;
                        current_delivery = reference.baseline;
                        outside = 0;
                    }
                }
            }
        }
        p = improved ? p0 : std::min(settings.max_relax, p * settings.relax_growth);
        operators.reward(which, move_score(outcome, settings));
        temperature *= settings.cooling;
    }
    return best;
}

template <typename Improve>
std::optional<Solution> multi_start(const Instance& instance, const Settings& settings,
                                    Deadline deadline, Stats& stats, Improve improve) {
    Random random(settings.seed);
    Operators construction(3, settings.construction_segment, settings.weight_reaction);
    std::optional<Solution> best;
    int no_improvement = 0;
    while (no_improvement < settings.max_no_improvement && Clock::now() < deadline) {
        ++stats.starts;
        const int which = construction.choose(random);
        ++stats.construction_uses[which];
        auto plan = construct(instance, which, settings, random, deadline);
        std::optional<Solution> solution;
        if (plan) {
            solution = improve(std::move(*plan), random);
        }
        const bool improved = solution && (!best || solution->timing.makespan() < best->timing.makespan() - tolerance);
        if (improved) {
            best = std::move(solution);
            no_improvement = 0;
        } else {
            ++no_improvement;
        }
        construction.reward(which, improved ? settings.construction_best_score : settings.construction_other_score);
    }
    return best;
}

}  // namespace

std::optional<Solution> search(GRBEnv& env, const Instance& instance, const Settings& settings,
                               Deadline deadline, Stats& stats) {
    return multi_start(instance, settings, deadline, stats, [&](Structure plan, Random& random)
                       -> std::optional<Solution> {
        const auto info = describe(instance, plan, true);
        if (!info || !prefix_feasible(instance, plan, *info)) {
            ++stats.screen_rejections;
            return std::nullopt;
        }
        auto solution = complete(env, instance, plan, {}, settings, random, deadline, stats);
        if (solution) {
            return local_search(env, instance, *solution, settings, random, deadline, stats);
        }
        return std::nullopt;
    });
}


}  // namespace fstsp
