#pragma once

#include <chrono>
#include <cstdint>
#include <iosfwd>
#include <optional>
#include <vector>

class GRBEnv;

namespace fstsp {

using Clock = std::chrono::steady_clock;
using Deadline = Clock::time_point;
constexpr double tolerance = 1e-7;

struct Customer {
    double x{}, y{}, weight{}, truck_service{}, drone_service{}, release{};
    bool delayed{}, rendezvous{};
};

struct Instance {
    std::vector<Customer> customers;
    std::vector<std::vector<double>> road, air;
    double delivery_capacity{}, resupply_capacity{}, delivery_endurance{};

    int end() const { return static_cast<int>(customers.size()) - 1; }
    int size() const { return end() - 1; }
    double flight_time(int launch, const std::vector<int>& visits, int recovery) const;
};

struct Sortie {
    int launch{}, recovery{};
    std::vector<int> visits;
};

struct Structure {
    std::vector<int> truck;
    std::vector<Sortie> sorties;
};

struct DeliveryInfo {
    std::vector<int> position, sortie_of, cutoff, rendezvous;
    std::vector<int> launch_pos, recovery_pos;
    std::vector<double> flight, payload;
};

struct Resupply {
    int node{};
    double departure{}, handoff{}, returns{};
    std::vector<int> parcels;
};

struct Schedule {
    std::vector<double> arrival, departure, launch, recovery;
    std::vector<Resupply> resupply;

    double makespan() const { return departure.back(); }
};

struct Solution {
    Structure delivery;
    Schedule timing;
    double baseline{};
};

struct Pair {
    int rendezvous{}, parcel{};
};

enum class FixedStatus { optimal, infeasible, time_limit };

struct FixedResult {
    FixedStatus status{FixedStatus::infeasible};
    std::optional<Solution> solution;
};

struct Stats {
    int starts{}, iterations{}, fixed_solves{}, screen_rejections{};
    int optimality_cuts{}, feasibility_cuts{}, warm_starts{}, warm_pairs{};
    int recovery_trials{}, screen_moves{}, forced_repairs{}, complete_improvements{};
    int construction_uses[3]{}, operator_uses[6]{};
};

enum class ClusterJoin { best_insertion, append };

struct Settings {
    std::uint64_t seed{};
    double seconds{};
    int max_no_improvement{};
    int max_infeasible_moves{};
    double cooling{};
    double noise{};
    double remove_fraction{};
    double max_relax{};
    double relax_growth{};
    int construction_segment{};
    int weight_segment{};
    double weight_reaction{};
    int max_clusters{};
    int cluster_iterations{};
    ClusterJoin cluster_join{};
    double shaw_exponent{};
    double initial_worsening{};
    double initial_acceptance{};
    double construction_best_score{};
    double construction_other_score{};
    double local_best_score{};
    double local_accepted_score{};
    double local_rejected_score{};
};

Instance read_instance(std::istream& input);
Structure read_structure(std::istream& input, const Instance& instance);
std::optional<DeliveryInfo> describe(const Instance& instance, const Structure& plan, bool complete);
bool prefix_feasible(const Instance& instance, const Structure& plan, const DeliveryInfo& info);
std::optional<Schedule> delivery_schedule(const Instance& instance, const Structure& plan,
                                          const DeliveryInfo& info, bool availability_bounds);
std::vector<Pair> assignments(const Solution& solution);
void verify_solution(const Instance& instance, const Solution& solution);
FixedResult solve_fixed(GRBEnv& env, const Instance& instance, const Structure& plan,
                       const std::vector<Pair>& warm, Deadline deadline, Stats& stats);
std::optional<Solution> search(GRBEnv& env, const Instance& instance, const Settings& settings,
                               Deadline deadline, Stats& stats);

}  // namespace fstsp
