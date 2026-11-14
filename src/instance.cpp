#include "fstsp.hpp"

#include <algorithm>
#include <cmath>
#include <istream>
#include <stdexcept>
#include <string>

namespace fstsp {
namespace {

bool less_equal(double a, double b) {
    return a <= b + tolerance;
}

void require(bool condition, const char* message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

bool nonnegative(double value) {
    return std::isfinite(value) && value >= 0;
}

}  // namespace

Instance read_instance(std::istream& input) {
    Instance result;
    int n, version;
    std::string magic;
    double truck_speed, drone_speed;
    require(static_cast<bool>(input >> magic >> version) && magic == "FSTSP_DR_LRP" && version == 1,
            "Expected FSTSP_DR_LRP 1.");
    require(static_cast<bool>(input >> n >> truck_speed >> drone_speed >> result.delivery_capacity >>
                             result.resupply_capacity >> result.delivery_endurance),
            "Expected n, speeds, capacities and delivery endurance.");
    require(n > 0 && nonnegative(truck_speed) && truck_speed > 0 &&
            nonnegative(drone_speed) && drone_speed > 0 &&
            nonnegative(result.delivery_capacity) && result.delivery_capacity > 0 &&
            nonnegative(result.resupply_capacity) && result.resupply_capacity > 0 &&
            nonnegative(result.delivery_endurance), "Invalid sizes, speeds, capacities or endurance.");
    result.customers.resize(static_cast<std::size_t>(n) + 2);
    auto& depot = result.customers.front();
    require(static_cast<bool>(input >> depot.x >> depot.y) &&
            std::isfinite(depot.x) && std::isfinite(depot.y), "Invalid depot coordinates.");
    for (int i = 1; i <= n; ++i) {
        auto& customer = result.customers[i];
        int id, delayed, rendezvous;
        require(static_cast<bool>(input >> id >> customer.x >> customer.y >> customer.weight >>
                                 customer.truck_service >> customer.drone_service >>
                                 customer.release >> delayed >> rendezvous), "Incomplete customer record.");
        require(id == i && std::isfinite(customer.x) && std::isfinite(customer.y) &&
                nonnegative(customer.weight) && customer.weight > 0 &&
                nonnegative(customer.truck_service) && nonnegative(customer.drone_service) &&
                nonnegative(customer.release) && (delayed == 0 || delayed == 1) &&
                (rendezvous == 0 || rendezvous == 1) && (delayed || customer.release == 0),
                "Invalid customer record (IDs must be consecutive, flags 0/1).");
        customer.delayed = delayed != 0;
        customer.rendezvous = rendezvous != 0;
    }
    result.customers.back() = depot;
    require(!(input >> magic), "Unexpected data after the customer records.");
    const auto dimension = result.customers.size();
    result.road.assign(dimension, std::vector<double>(dimension));
    result.air = result.road;
    for (int i = 0; i <= result.end(); ++i) {
        for (int j = i + 1; j <= result.end(); ++j) {
            const double dx = result.customers[i].x - result.customers[j].x;
            const double dy = result.customers[i].y - result.customers[j].y;
            const double road = (std::abs(dx) + std::abs(dy)) / truck_speed;
            const double air = std::hypot(dx, dy) / drone_speed;
            require(std::isfinite(road) && std::isfinite(air), "Non-finite travel time.");
            require(less_equal(air, road), "Drone travel time must not exceed truck travel time.");
            result.road[i][j] = result.road[j][i] = road;
            result.air[i][j] = result.air[j][i] = air;
        }
    }
    return result;
}

Structure read_structure(std::istream& input, const Instance& instance) {
    Structure result;
    std::string tag;
    int count;
    require(static_cast<bool>(input >> tag >> count) && tag == "truck" &&
            count >= 3 && count <= instance.end() + 1, "Expected truck <count> <nodes>.");
    result.truck.resize(count);
    for (int& node : result.truck) {
        require(static_cast<bool>(input >> node), "Missing truck node.");
    }
    require(static_cast<bool>(input >> tag >> count) && tag == "sorties" &&
            count >= 0 && count <= instance.size(), "Expected sorties <count>.");
    result.sorties.resize(count);
    for (auto& sortie : result.sorties) {
        int visits;
        require(static_cast<bool>(input >> sortie.launch >> sortie.recovery >> visits) &&
                visits >= 1 && visits <= instance.size(), "Expected launch recovery <count> <visits>.");
        sortie.visits.resize(visits);
        for (int& node : sortie.visits) {
            require(static_cast<bool>(input >> node), "Missing drone customer.");
        }
    }
    require(!(input >> tag), "Unexpected data after the sorties.");
    return result;
}

double Instance::flight_time(int launch, const std::vector<int>& visits, int recovery) const {
    double duration = 0;
    int previous = launch;
    for (const int node : visits) {
        duration += air[previous][node] + customers[node].drone_service;
        previous = node;
    }
    return duration + air[previous][recovery];
}

std::optional<DeliveryInfo> describe(const Instance& instance, const Structure& plan, bool complete) {
    const int end = instance.end();
    if (plan.truck.size() < (complete ? 3u : 2u) ||
        plan.truck.front() != 0 || plan.truck.back() != end) {
        return std::nullopt;
    }
    DeliveryInfo result;
    result.position.assign(instance.customers.size(), -1);
    result.sortie_of = result.position;
    result.cutoff = result.position;
    for (int k = 0; k < static_cast<int>(plan.truck.size()); ++k) {
        const int node = plan.truck[k];
        if (node < 0 || node > end || result.position[node] >= 0) {
            return std::nullopt;
        }
        result.position[node] = result.cutoff[node] = k;
        if (instance.customers[node].rendezvous) {
            result.rendezvous.push_back(k);
        }
    }
    int last_recovery = 0;
    for (int s = 0; s < static_cast<int>(plan.sorties.size()); ++s) {
        const auto& sortie = plan.sorties[s];
        if (sortie.launch < 0 || sortie.launch > end || sortie.recovery < 0 ||
            sortie.recovery > end || sortie.visits.empty()) {
            return std::nullopt;
        }
        const int launch = result.position[sortie.launch];
        const int recovery = result.position[sortie.recovery];
        if (launch < last_recovery || recovery <= launch) {
            return std::nullopt;
        }
        last_recovery = recovery;
        double load = 0;
        for (const int node : sortie.visits) {
            if (node <= 0 || node >= end || result.cutoff[node] >= 0) {
                return std::nullopt;
            }
            load += instance.customers[node].weight;
            result.sortie_of[node] = s;
            result.cutoff[node] = launch;
        }
        const double flight = instance.flight_time(sortie.launch, sortie.visits, sortie.recovery);
        if (!less_equal(load, instance.delivery_capacity) || !less_equal(flight, instance.delivery_endurance)) {
            return std::nullopt;
        }
        result.launch_pos.push_back(launch);
        result.recovery_pos.push_back(recovery);
        result.flight.push_back(flight);
        result.payload.push_back(load);
    }
    if (complete && std::find(result.cutoff.begin(), result.cutoff.end(), -1) != result.cutoff.end()) {
        return std::nullopt;
    }
    return result;
}

bool prefix_feasible(const Instance& instance, const Structure& plan, const DeliveryInfo& info) {
    // (47) in its equivalent unrounded form.
    std::vector<double> due(plan.truck.size());
    for (int j = 1; j <= instance.size(); ++j) {
        if (instance.customers[j].delayed && info.cutoff[j] >= 0) {
            due[info.cutoff[j]] += instance.customers[j].weight;
        }
    }
    int stops = 0;
    double demand = 0;
    for (int k = 0; k + 1 < static_cast<int>(plan.truck.size()); ++k) {
        if (instance.customers[plan.truck[k]].rendezvous) {
            ++stops;
        }
        demand += due[k];
        if (!less_equal(demand, stops * instance.resupply_capacity)) {
            return false;
        }
    }
    return true;
}

std::vector<Pair> assignments(const Solution& solution) {
    std::vector<Pair> result;
    for (const auto& trip : solution.timing.resupply) {
        for (const int parcel : trip.parcels) {
            result.push_back({trip.node, parcel});
        }
    }
    return result;
}

void verify_solution(const Instance& instance, const Solution& solution) {
    const auto& plan = solution.delivery;
    const auto& t = solution.timing;
    const auto info = describe(instance, plan, true);
    require(info.has_value(), "Output delivery structure is infeasible.");
    require(t.arrival.size() == plan.truck.size() && t.departure.size() == plan.truck.size() &&
            t.launch.size() == plan.sorties.size() && t.recovery.size() == plan.sorties.size(),
            "Output timetable size does not match the delivery structure.");
    require(std::abs(t.arrival.front()) <= tolerance && less_equal(0, t.departure.front()),
            "Invalid depot departure.");
    for (int k = 1; k < static_cast<int>(plan.truck.size()); ++k) {
        const int i = plan.truck[k - 1], j = plan.truck[k];
        require(std::abs(t.arrival[k] - t.departure[k - 1] - instance.road[i][j]) <= tolerance &&
                less_equal(t.arrival[k] + instance.customers[j].truck_service, t.departure[k]),
                "Invalid truck travel or service time.");
    }
    for (int s = 0; s < static_cast<int>(plan.sorties.size()); ++s) {
        const int l = info->launch_pos[s], b = info->recovery_pos[s];
        const double duration = t.recovery[s] - t.launch[s];
        require(less_equal(t.arrival[l], t.launch[s]) && less_equal(t.launch[s], t.departure[l]) &&
                less_equal(t.arrival[b], t.recovery[s]) && less_equal(t.recovery[s], t.departure[b]) &&
                less_equal(info->flight[s], duration) && less_equal(duration, instance.delivery_endurance) &&
                (s == 0 || less_equal(t.recovery[s - 1], t.launch[s])),
                "Invalid delivery launch, recovery or endurance.");
    }
    std::vector<int> supplied(instance.customers.size());
    double last_return = 0;
    int last_position = 0;
    for (const auto& trip : t.resupply) {
        require(trip.node > 0 && trip.node < instance.end() &&
                instance.customers[trip.node].rendezvous && !trip.parcels.empty(), "Invalid resupply stop.");
        const int position = info->position[trip.node];
        const double tau = instance.air[0][trip.node];
        require(position > last_position && less_equal(last_return, trip.departure) &&
                std::abs(trip.handoff - trip.departure - tau) <= tolerance &&
                std::abs(trip.returns - trip.handoff - tau) <= tolerance &&
                less_equal(t.arrival[position], trip.handoff) && less_equal(trip.handoff, t.departure[position]),
                "Invalid resupply order or handoff window.");
        double load = 0;
        for (const int j : trip.parcels) {
            require(j > 0 && j < instance.end() && instance.customers[j].delayed &&
                    ++supplied[j] == 1 && position <= info->cutoff[j] &&
                    less_equal(instance.customers[j].release, trip.departure),
                    "Invalid resupply parcel, release or handoff position.");
            const int s = info->sortie_of[j];
            if (s >= 0) {
                require(less_equal(trip.handoff, t.launch[s]), "Delivery launched before handoff.");
            }
            if (j == trip.node) {
                require(less_equal(trip.handoff + instance.customers[j].truck_service, t.departure[position]),
                        "Truck service starts before the same-stop handoff.");
            }
            load += instance.customers[j].weight;
        }
        require(less_equal(load, instance.resupply_capacity), "Resupply payload exceeds capacity.");
        last_return = trip.returns;
        last_position = position;
    }
    for (int j = 1; j <= instance.size(); ++j) {
        require(supplied[j] == static_cast<int>(instance.customers[j].delayed), "Missing delayed parcel.");
    }
}

}  // namespace fstsp
