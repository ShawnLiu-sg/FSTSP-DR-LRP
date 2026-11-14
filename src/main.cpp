#include "fstsp.hpp"

#include "gurobi_c++.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <set>
#include <utility>
#include <vector>

namespace {

using namespace fstsp;

struct Command {
    Settings settings;
    std::string instance, fixed;
};

void usage() {
    std::cout << "Usage: fstsp_dr_lrp INSTANCE --config FILE [--fixed STRUCTURE] [--seed N] [--seconds S]\n"
                 "Fill in algorithm.params locally. Command-line switches override file values.\n";
}

double real(const std::string& value) {
    std::size_t end;
    const double result = std::stod(value, &end);
    if (end != value.size() || !std::isfinite(result)) {
        throw std::runtime_error("Invalid numeric option: " + value);
    }
    return result;
}

std::uint64_t natural(const std::string& value) {
    if (value.empty() || value.find_first_not_of("0123456789") != std::string::npos) {
        throw std::runtime_error("Expected a nonnegative integer: " + value);
    }
    return std::stoull(value);
}

int positive_integer(const std::string& value) {
    const auto number = natural(value);
    if (number == 0 || number > static_cast<std::uint64_t>(std::numeric_limits<int>::max())) {
        throw std::runtime_error("Expected a positive integer: " + value);
    }
    return static_cast<int>(number);
}

void apply_option(Command& command, const std::string& option, const std::string& value) {
    auto& s = command.settings;
    if (option == "--fixed") {
        command.fixed = value;
    } else if (option == "--seed") {
        s.seed = natural(value);
    } else if (option == "--seconds") {
        s.seconds = real(value);
        if (s.seconds <= 0) {
            throw std::runtime_error("--seconds must be positive.");
        }
    } else if (option == "--no-update") {
        s.max_no_improvement = positive_integer(value);
    } else if (option == "--max-infeasible") {
        s.max_infeasible_moves = positive_integer(value);
    } else if (option == "--noise") {
        s.noise = real(value);
    } else if (option == "--remove") {
        s.remove_fraction = real(value);
    } else if (option == "--max-relax") {
        s.max_relax = real(value);
    } else if (option == "--relax-growth") {
        s.relax_growth = real(value);
    } else if (option == "--cooling") {
        s.cooling = real(value);
    } else if (option == "--construction-segment") {
        s.construction_segment = positive_integer(value);
    } else if (option == "--weight-segment") {
        s.weight_segment = positive_integer(value);
    } else if (option == "--weight-reaction") {
        s.weight_reaction = real(value);
    } else if (option == "--max-clusters") {
        s.max_clusters = positive_integer(value);
    } else if (option == "--cluster-iterations") {
        s.cluster_iterations = positive_integer(value);
    } else if (option == "--cluster-join") {
        if (value == "best-insertion") {
            s.cluster_join = ClusterJoin::best_insertion;
        } else if (value == "append") {
            s.cluster_join = ClusterJoin::append;
        } else {
            throw std::runtime_error("--cluster-join expects best-insertion or append.");
        }
    } else if (option == "--shaw-exponent") {
        s.shaw_exponent = real(value);
    } else if (option == "--initial-worsening") {
        s.initial_worsening = real(value);
    } else if (option == "--initial-acceptance") {
        s.initial_acceptance = real(value);
    } else if (option == "--construction-best-score") {
        s.construction_best_score = real(value);
    } else if (option == "--construction-other-score") {
        s.construction_other_score = real(value);
    } else if (option == "--local-best-score") {
        s.local_best_score = real(value);
    } else if (option == "--local-accepted-score") {
        s.local_accepted_score = real(value);
    } else if (option == "--local-rejected-score") {
        s.local_rejected_score = real(value);
    } else {
        throw std::runtime_error("Unknown option: " + option);
    }
}

void read_config(Command& command, const std::string& path, std::set<std::string>& supplied) {
    std::ifstream file(path);
    if (!file) {
        throw std::runtime_error("Cannot read parameter file: " + path);
    }
    std::string line;
    while (std::getline(file, line)) {
        if (const auto comment = line.find('#'); comment != std::string::npos) {
            line.erase(comment);
        }
        std::replace(line.begin(), line.end(), '=', ' ');
        std::istringstream row(line);
        std::string key, value, extra;
        if (!(row >> key)) {
            continue;
        }
        // A blank template row leaves the option to be supplied by the caller.
        if (!(row >> value)) {
            continue;
        }
        if (row >> extra || key == "fixed") {
            throw std::runtime_error("Expected parameter key=value in " + path + ": " + line);
        }
        apply_option(command, "--" + key, value);
        supplied.insert(key);
    }
}

Command parse(int argc, char** argv) {
    Command command;
    command.instance = argv[1];
    std::vector<std::pair<std::string, std::string>> overrides;
    std::set<std::string> supplied;
    for (int k = 2; k < argc; k += 2) {
        const std::string option = argv[k];
        if (k + 1 == argc) {
            throw std::runtime_error("Missing value for " + option);
        }
        if (option == "--config") {
            read_config(command, argv[k + 1], supplied);
        } else {
            overrides.emplace_back(option, argv[k + 1]);
        }
    }
    for (const auto& [option, value] : overrides) {
        apply_option(command, option, value);
        supplied.insert(option.substr(2));
    }
    for (const char* option : {"seed", "seconds", "no-update", "max-infeasible",
                               "noise", "remove", "max-relax", "relax-growth",
                               "cooling", "construction-segment", "weight-segment", "weight-reaction",
                               "max-clusters", "cluster-iterations", "cluster-join", "shaw-exponent",
                               "initial-worsening", "initial-acceptance", "construction-best-score", "construction-other-score",
                               "local-best-score", "local-accepted-score", "local-rejected-score"}) {
        if (supplied.count(option) == 0 && (command.fixed.empty() || std::string(option) == "seconds")) {
            throw std::runtime_error(std::string("Supply a value for ") + option + " in the parameter file or on the command line.");
        }
    }
    auto& s = command.settings;
    if (command.fixed.empty() && (s.noise < 0 || s.noise > 1 || s.remove_fraction <= 0 || s.remove_fraction > 1 ||
        s.max_relax < 0 || s.max_relax > 1 || s.relax_growth <= 1 || s.cooling <= 0 || s.cooling > 1 ||
        s.weight_reaction < 0 || s.weight_reaction > 1 || s.shaw_exponent <= 0 ||
        s.initial_worsening <= 0 || s.initial_acceptance <= 0 || s.initial_acceptance >= 1 ||
        s.construction_best_score <= 0 || s.construction_other_score <= 0 ||
        s.local_best_score <= 0 || s.local_accepted_score <= 0 || s.local_rejected_score <= 0)) {
        throw std::runtime_error("Search parameters are outside their allowed ranges.");
    }
    return command;
}

void print_nodes(const std::vector<int>& nodes) {
    std::cout << '[';
    for (std::size_t k = 0; k < nodes.size(); ++k) {
        std::cout << (k ? "," : "") << nodes[k];
    }
    std::cout << ']';
}

void print_delivery(const Instance& instance, const Solution& solution, double start, const char* indent) {
    const auto& x = solution.delivery;
    const auto& t = solution.timing;
    std::cout << ",\n" << indent << "\"truck\": [";
    for (std::size_t k = 0; k < x.truck.size(); ++k) {
        const int node = x.truck[k];
        std::cout << (k ? "," : "") << '\n' << indent << "  {\"node\":" << node
                  << ",\"arrival\":" << start + t.arrival[k]
                  << ",\"service_start\":" << start + (t.departure[k] - instance.customers[node].truck_service)
                  << ",\"departure\":" << start + t.departure[k] << '}';
    }
    std::cout << '\n' << indent << "],\n" << indent << "\"delivery_sorties\": [";
    for (std::size_t sidx = 0; sidx < x.sorties.size(); ++sidx) {
        const auto& sortie = x.sorties[sidx];
        std::cout << (sidx ? "," : "") << '\n' << indent << "  {\"launch\":" << sortie.launch
                  << ",\"takeoff\":" << start + t.launch[sidx] << ",\"visits\":";
        print_nodes(sortie.visits);
        std::cout << ",\"services\":[";
        int last = sortie.launch;
        double next = t.launch[sidx];
        for (std::size_t j = 0; j < sortie.visits.size(); ++j) {
            const int node = sortie.visits[j];
            const double arrive = next + instance.air[last][node];
            next = arrive + instance.customers[node].drone_service;
            std::cout << (j ? "," : "") << "{\"node\":" << node << ",\"arrival\":" << start + arrive
                      << ",\"service_end\":" << start + next << '}';
            last = node;
        }
        std::cout << "],\"recovery\":" << sortie.recovery << ",\"land\":" << start + t.recovery[sidx]
                  << ",\"endurance_used\":" << t.recovery[sidx] - t.launch[sidx] << '}';
    }
    std::cout << '\n' << indent << ']';
}

void print_header(const char* status, const Settings& settings, double elapsed) {
    std::cout << std::setprecision(15)
              << "{\n  \"status\": \"" << status << "\","
              << "\n  \"seed\": " << settings.seed
              << ",\n  \"time_limit_seconds\": " << settings.seconds
              << ",\n  \"elapsed_seconds\": " << elapsed
              << ",\n  \"algorithm_parameters\": {"
              << "\"no-update\":" << settings.max_no_improvement
              << ",\"max-infeasible\":" << settings.max_infeasible_moves
              << ",\"noise\":" << settings.noise
              << ",\"remove\":" << settings.remove_fraction
              << ",\"max-relax\":" << settings.max_relax
              << ",\"relax-growth\":" << settings.relax_growth
              << ",\"cooling\":" << settings.cooling
              << ",\"construction-segment\":" << settings.construction_segment
              << ",\"weight-segment\":" << settings.weight_segment
              << ",\"weight-reaction\":" << settings.weight_reaction
              << ",\"max-clusters\":" << settings.max_clusters
              << ",\"cluster-iterations\":" << settings.cluster_iterations
              << ",\"cluster-join\":\""
              << (settings.cluster_join == ClusterJoin::append ? "append" : "best-insertion") << '"'
              << ",\"shaw-exponent\":" << settings.shaw_exponent
              << ",\"initial-worsening\":" << settings.initial_worsening
              << ",\"initial-acceptance\":" << settings.initial_acceptance
              << ",\"construction-best-score\":" << settings.construction_best_score
              << ",\"construction-other-score\":" << settings.construction_other_score
              << ",\"local-best-score\":" << settings.local_best_score
              << ",\"local-accepted-score\":" << settings.local_accepted_score
              << ",\"local-rejected-score\":" << settings.local_rejected_score << '}';
}

void print_statistics(const Stats& stats) {
    std::cout << ",\n  \"statistics\": {\"starts\":" << stats.starts << ",\"iterations\":" << stats.iterations
              << ",\"fixed_solves\":" << stats.fixed_solves
              << ",\"screen_rejections\":" << stats.screen_rejections
              << ",\"optimality_cuts\":" << stats.optimality_cuts << ",\"feasibility_cuts\":" << stats.feasibility_cuts
              << ",\"warm_starts\":" << stats.warm_starts << ",\"warm_pairs\":" << stats.warm_pairs
              << ",\"recovery_trials\":" << stats.recovery_trials << ",\"screen_moves\":" << stats.screen_moves
              << ",\"forced_repairs\":" << stats.forced_repairs << ",\"complete_improvements\":"
              << stats.complete_improvements << ",\"construction_uses\":";
    print_nodes(std::vector<int>(std::begin(stats.construction_uses), std::end(stats.construction_uses)));
    std::cout << ",\"operator_uses\":";
    print_nodes(std::vector<int>(std::begin(stats.operator_uses), std::end(stats.operator_uses)));
    std::cout << "}\n}\n";
}

void print_result(const Instance& instance, const std::optional<Solution>& solution, const char* status,
                  const Settings& s, const Stats& stats, double elapsed) {
    print_header(status, s, elapsed);
    if (solution) {
        const auto& t = solution->timing;
        std::cout << ",\n  \"makespan\": " << t.makespan()
                  << ",\n  \"delivery_only_baseline\": " << solution->baseline
                  << ",\n  \"resupply_delay\": " << t.makespan() - solution->baseline;
        print_delivery(instance, *solution, 0, "  ");
        std::cout << ",\n  \"resupply_sorties\": [";
        for (std::size_t i = 0; i < t.resupply.size(); ++i) {
            const auto& trip = t.resupply[i];
            std::cout << (i ? "," : "") << "\n    {\"node\":" << trip.node << ",\"takeoff\":" << trip.departure
                      << ",\"handoff\":" << trip.handoff << ",\"return\":" << trip.returns << ",\"parcels\":";
            print_nodes(trip.parcels);
            std::cout << '}';
        }
        std::cout << "\n  ]";
    }
    print_statistics(stats);
}

}  // namespace

int main(int argc, char** argv) {
    if (argc == 1 || std::string(argv[1]) == "--help") {
        usage();
        return argc == 1 ? 1 : 0;
    }
    try {
        auto command = parse(argc, argv);
        std::ifstream input(command.instance);
        if (!input) {
            throw std::runtime_error("Cannot read instance: " + command.instance);
        }
        const auto instance = read_instance(input);
        const auto budget = std::chrono::duration_cast<Clock::duration>(
            std::chrono::duration<double>(command.settings.seconds));
        GRBEnv env(true);
        env.set(GRB_IntParam_OutputFlag, 0);
        env.start();
        const auto begin = Clock::now();
        const auto deadline = begin + budget;
        Stats stats;
        std::optional<Solution> solution;
        const char* status = "no_feasible_solution_found";
        if (command.fixed.empty()) {
            solution = search(env, instance, command.settings, deadline, stats);
            status = solution ? "feasible" : "no_feasible_solution_found";
        } else {
            std::ifstream structure(command.fixed);
            if (!structure) {
                throw std::runtime_error("Cannot read structure: " + command.fixed);
            }
            auto result = solve_fixed(env, instance, read_structure(structure, instance), {}, deadline, stats);
            solution = std::move(result.solution);
            switch (result.status) {
                case FixedStatus::optimal: status = "optimal"; break;
                case FixedStatus::infeasible: status = "infeasible"; break;
                case FixedStatus::time_limit: status = "time_limit"; break;
            }
        }
        if (solution) {
            verify_solution(instance, *solution);
        }
        print_result(instance, solution, status, command.settings, stats,
                     std::chrono::duration<double>(Clock::now() - begin).count());
        return 0;
    } catch (const GRBException& error) {
        std::cerr << "Gurobi error " << error.getErrorCode() << ": " << error.getMessage() << '\n';
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
    }
    return 1;
}
