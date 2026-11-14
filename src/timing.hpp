#pragma once

#include "fstsp.hpp"

#include <vector>

namespace fstsp {

// An arc from u to v represents t[v] - t[u] >= rhs.
// A term with complement=true multiplies (1 - x) instead of x.
struct Term {
    int variable{};
    double coefficient{};
    bool complement{};
};

struct TimeArc {
    int from{}, to{};
    double base{};
    std::vector<Term> terms;
};

struct TimeAnswer {
    bool feasible{};
    std::vector<double> times;
    std::vector<int> certificate;
};

struct TimeGraph {
    explicit TimeGraph(int size, int sorties, int stops = 0);

    int arrival(int k) const { return 1 + 2 * k; }
    int departure(int k) const { return arrival(k) + 1; }
    int launch(int s) const { return 1 + 2 * (truck_count + s); }
    int recovery(int s) const { return launch(s) + 1; }
    int resupply(int i) const { return 1 + 2 * (truck_count + sortie_count) + i; }
    int finish() const { return departure(truck_count - 1); }

    void add(int from, int to, double rhs, std::vector<Term> terms = {});
    void delivery_rows(const Instance& instance, const Structure& plan, const DeliveryInfo& info,
                       bool availability_bounds);
    TimeAnswer solve(const std::vector<int>& binary, bool get_certificate = false) const;
    Schedule timetable(const std::vector<double>& times) const;

    int truck_count, sortie_count, vertices;
    std::vector<TimeArc> arcs;
};

}  // namespace fstsp
