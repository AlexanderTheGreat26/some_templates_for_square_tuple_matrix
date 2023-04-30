#include <iostream>
#include <tuple>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <iterator>
#include <memory>
#include <algorithm>
#include <string>
#include <array>


const double max_cation_length ();



typedef std::tuple<double, double, double> point;
typedef std::vector<point> frame;
typedef std::vector<frame> frames;


point geometric_center_between_O (point & O_1, point & O_2);


int main() {



    // First of all we have to find geometric centers of cations. Note that we have Eigen in the system.


    return 0;
}


template<size_t Is = 0, typename... Tp>
void vector_scalar_multiplication (std::tuple<Tp...>& vector, const double & lambda, std::tuple<Tp...>& result) {
    std::get<Is>(result) = std::get<Is>(vector) * lambda; // May be better to make a recursion, but this more universal.
    if constexpr(Is + 1 != sizeof...(Tp))
        vector_scalar_multiplication<Is + 1>(vector, lambda, result);
}


template<size_t Is = 0, typename... Tp>
void vector_offset (std::tuple<Tp...>& vector, std::tuple<Tp...>& frame_of_reference, std::tuple<Tp...>& result) {
    std::get<Is>(result) = std::get<Is>(vector) + std::get<Is>(frame_of_reference);
    if constexpr(Is + 1 != sizeof...(Tp))
        vector_offset<Is + 1>(vector, frame_of_reference, result);
}


point mass_center (std::vector<point> & coordinates, std::vector<double> & masses) {
    point R, tmp, mR = std::make_tuple(0, 0, 0);
    double M = 0;
    for (int i = 0; i < coordinates.size(); ++i) {
        vector_scalar_multiplication(coordinates[i], masses[i], tmp);
        vector_offset(tmp, mR, mR); // mR = tmp + mR
        M += masses[i];
    }
    vector_scalar_multiplication(mR, 1.0/M, R);
    return R;
}


point geometric_center_between_O (point & O_1, point & O_2) {
    point result;
    vector_offset(O_1, O_2, result);
    vector_scalar_multiplication(result, 0.5, result);
    return result;
}


template<typename T, size_t... Is>
double distance_impl (T const& t, T const& t1, std::index_sequence<Is...>, std::index_sequence<Is...>) {
    return (std::sqrt((std::pow(std::get<Is>(t) - std::get<Is>(t1), 2) + ...)));
}

template <class Tuple>
double distance (const Tuple & t, const Tuple & t1) {
    constexpr auto size = std::tuple_size<Tuple>{};
    return distance_impl(t, t1, std::make_index_sequence<size>{}, std::make_index_sequence<size>{});
}


std::vector<double> distances_in_frame (frame & atoms) {
    std::vector<double> results;
    for (int i = 0; i < atoms.size(); ++i)
        for (int j = 0; j < i; ++j)
            results.emplace_back(distance(atoms[i], atoms[j]));
    return results;
}


std::vector<point> potential_cations (frame & oxygens, const double & max_distance) {
    std::vector<point> result;
    for (int i = 0; i < oxygens.size(); ++i) {
        for (int j = 0; j < i; ++j)
            if (distance(oxygens[i], oxygens[j]) <= max_distance) {
                result.emplace_back(oxygens[i]);
                result.emplace_back(oxygens[j]);
            }
    }
    return result;
}


bool pair_of_interest (point & O_1, point & O_2, frame & free_proton, ) {

}


