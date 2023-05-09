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


const size_t dim = std::tuple_size<point>{};


point geometric_center_between_O (point & O_1, point & O_2);


int main() {



    // First of all we have to find geometric centers of cations. Note that we have Eigen in the system.


    return 0;
}


template<size_t Is = 0, typename... Tp>
void vector_scalar_multiplication (const std::tuple<Tp...>& vector, const double & lambda, std::tuple<Tp...>& result) {
    std::get<Is>(result) = std::get<Is>(vector) * lambda; // May be better to make a recursion, but this more universal.
    if constexpr(Is + 1 != sizeof...(Tp))
        vector_scalar_multiplication<Is + 1>(vector, lambda, result);
}


template<size_t Is = 0, typename... Tp>
void vector_offset (const std::tuple<Tp...>& vector, std::tuple<Tp...>& frame_of_reference, std::tuple<Tp...>& result) {
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
    for (int i = 0; i < oxygens.size(); ++i)
        for (int j = 0; j < i; ++j)
            if (distance(oxygens[i], oxygens[j]) <= max_distance) {
                result.emplace_back(oxygens[i]);
                result.emplace_back(oxygens[j]);
            }
    return result;
}  // Запусти покадрово и найдёшь максимальный размер катиона,


bool is_Zundel (point & O_1, point & O_2, frame & free_protons, const double & deviation) {
    bool ans = false;
    point center = geometric_center_between_O(O_2, O_1);
    for (const auto & free_proton : free_protons)
        if (distance(free_proton, center) <= deviation) {
            ans = true;
            break;
        }
    return ans;
}


template<typename T, size_t... Is>
double scalar_vector_multiplication_impl (T const& a, T const& b, std::index_sequence<Is...>, std::index_sequence<Is...>) {
    return ((std::get<Is>(a) * std::get<Is>(b)) + ...);
}

template <class Tuple>
double scalar_vector_multiplication (const Tuple & a, const Tuple & b) {
    constexpr auto size = std::tuple_size<Tuple>{};
    return scalar_vector_multiplication_impl(a, b, std::make_index_sequence<size>{}, std::make_index_sequence<size>{});
}


template <class Tuple>
double vector_abs (const Tuple & t) {
    constexpr auto size = std::tuple_size<Tuple>{};
    return std::sqrt(scalar_vector_multiplication_impl(t, t, std::make_index_sequence<size>{}, std::make_index_sequence<size>{}));
}


template <class Tuple>
double cos_angle_between_vectors  (const Tuple & a, const Tuple & b) {
    return (scalar_vector_multiplication(a, b) / vector_abs(a) / vector_abs(b));
}


point vector_vector_multiplication (const point & a, const point & b) {
    double x = std::get<1>(a)*std::get<2>(b) - std::get<2>(a)*std::get<1>(b);
    double y = std::get<0>(a)*std::get<2>(b) - std::get<2>(a)*std::get<0>(b);
    double z = std::get<0>(a)*std::get<1>(b) - std::get<1>(a)*std::get<0>(b);
    return std::make_tuple(x, -y, z);
}


double sin_angle_between_vectors  (const point & a, const point & b) {
    return vector_abs(vector_vector_multiplication(a, b)) / vector_abs(a) / vector_abs(b);
}








std::vector<point> frame_of_reference_rotation (point & direction_vector, const std::vector<point> & basis_set) {
    std::vector<point> new_basis_set;

//    double cos_a, cos_b, cos_g, sin_a, sin_b, sin_g;
//    cos_a = cos_angle_between_vectors(direction_vector, basis_set[0]);
//    cos_b = cos_angle_between_vectors(direction_vector, basis_set[1]);
//    cos_g = cos_angle_between_vectors(direction_vector, basis_set[2]);
//    sin_a = sin_angle_between_vectors(direction_vector, basis_set[0]);
//    sin_b = sin_angle_between_vectors(direction_vector, basis_set[1]);
//    sin_g = sin_angle_between_vectors(direction_vector, basis_set[2]);

    std::vector<point> M; // Rotation matrix
    M.emplace_back(cos_b*cos_g,                     -sin_a*cos_b,                      sin_b);
    M.emplace_back(sin_a*sin_b*cos_g + sin_g*cos_a, -sin_a*sin_b*sin_g + cos_a*cos_g, -sin_a*cos_b);
    M.emplace_back(sin_a*sin_g - sin_b*cos_a*cos_g,  sin_a*cos_g + sin_b*sin_g*cos_a,  cos_a*cos_b);

    new_basis_set.reserve(basis_set.size());
    for (const auto & vector : basis_set)
        new_basis_set.emplace_back(scalar_vector_multiplication(M[0], vector),
                                   scalar_vector_multiplication(M[1], vector),
                                   scalar_vector_multiplication(M[2], vector)); // Чёт говной воняет.

    return new_basis_set;
}


bool is_cation (point & O_1, point & O_2, frame & free_protons, const double & deviation) {
    bool ans = false;
    point ref, O_2_in_O1_ref;
    vector_scalar_multiplication(O_1, -1.0, ref);
    vector_offset(O_2, ref, O_2_in_O1_ref);
    const std::vector<point> default_basis_set = {std::make_tuple(1, 0, 0),
                                                  std::make_tuple(0, 1, 0),
                                                  std::make_tuple(0, 0, 1)};

    std::vector<point> new_basis_set = std::move(frame_of_reference_rotation(O_2_in_O1_ref, default_basis_set));

    for (int i = 0; i < free_protons.size(); ++i) {
        vector_offset(, ref, O_2_in_O1_ref);
    }

    return ans;
}


bool pair_of_interest (point & O_1, point & O_2, frame & free_proton, const double & deviation, const std::string & type) {
    if (type == "Zundel") {
        point geometric_center = geometric_center_between_O(O_1,O_2);
    }

}


point direction_of_rotation (point & vector, std::vector<point> & basis_set) {
   point result;

    // We need to create a projection to xOy of our vector
    point xOy_projection = std::make_tuple(std::get<0>(vector), std::get<1>(vector), 0);

    double cos_a, cos_b, cos_g, sin_a, sin_b, sin_g;

    cos_a = cos_angle_between_vectors(xOy_projection, basis_set[0]);
    sin_a = sin_angle_between_vectors(xOy_projection, basis_set[0]);

    cos_b = cos_angle_between_vectors(xOy_projection, basis_set[1]);
    sin_b = sin_angle_between_vectors(xOy_projection, basis_set[1]);

    cos_g = cos_angle_between_vectors(vector, basis_set[2]);
    sin_g = sin_angle_between_vectors(vector, basis_set[2]);
    return result;
}


template <class Tuple>
auto linear_vector_sum  (const Tuple & a, const double & lambda_1, const Tuple & b, const double & lambda_2) {
    auto buf_a = a;
    auto buf_b = b;
    vector_scalar_multiplication(a, lambda_1, buf_a);
    vector_scalar_multiplication(b, lambda_2, buf_b);
    vector_offset(buf_b, buf_a, buf_b);
    return buf_b;
}


template<size_t Is = 0, typename... Tp>
void forward (std::vector<std::tuple<Tp...>> & invertible_matrix, std::vector<std::tuple<Tp...>> & identity_matrix) {
    for (int i = Is+1; i < invertible_matrix.size(); ++i) {
        double lambda = -std::get<Is>(invertible_matrix[i]) / std::get<Is>(invertible_matrix[Is]);
        invertible_matrix[i] = linear_vector_sum(invertible_matrix[Is], lambda, invertible_matrix[i], 1.0);
        identity_matrix[i] = linear_vector_sum(identity_matrix[Is], lambda, identity_matrix[i], 1.0);
    }
    if constexpr (Is + 1 != sizeof...(Tp))
        forward<Is + 1>(invertible_matrix, identity_matrix);
}


template<size_t Is = 0, typename... Tp, size_t N = sizeof...(Tp)-1>
void back (std::vector<std::tuple<Tp...>> & invertible_matrix, std::vector<std::tuple<Tp...>> & identity_matrix) {
    for (int i = N-Is; i > 0; --i) {
        double lambda = -std::get<N-Is>(invertible_matrix[i-1]) / std::get<N-Is>(invertible_matrix[N-Is]);
        invertible_matrix[i-1] = linear_vector_sum(invertible_matrix[N-Is], lambda, invertible_matrix[i-1], 1.0);
        identity_matrix[i-1] = linear_vector_sum(identity_matrix[N-Is], lambda, identity_matrix[i-1], 1.0);
    }
    if constexpr (N-Is > 0)
        back<Is + 1>(invertible_matrix, identity_matrix);
}


template<size_t Is = 0, typename... Tp>
void final_step (std::vector<std::tuple<Tp...>> & invertible_matrix, std::vector<std::tuple<Tp...>> & identity_matrix) {
    vector_scalar_multiplication(identity_matrix[Is], 1.0 / std::get<Is>(invertible_matrix[Is]), identity_matrix[Is]);
    if constexpr (Is + 1 != sizeof...(Tp))
        final_step<Is + 1>(invertible_matrix, identity_matrix);
}


template<size_t Is = 0, typename... Tp>
void identity_matrix (std::vector<std::tuple<Tp...>> & E) {
    for (int i = 0; i < E.size(); ++i)
        std::get<Is>(E[i]) = (i == Is) ? 1 : 0;
    if constexpr (Is + 1 != sizeof...(Tp))
        identity_matrix<Is + 1>(E);
}


template<size_t Is = 0, typename... Tp>
auto Gaussian_elimination (std::vector<std::tuple<Tp...>> invertible_matrix) {
    auto E = invertible_matrix;
    identity_matrix(E);
    forward(invertible_matrix, E);
    back(invertible_matrix, E);
    final_step(invertible_matrix, E);
    return E;
}


template <typename F, size_t... Is>
auto gen_tuple_impl(F & func, std::index_sequence<Is...> ) {
    return std::make_tuple(func(Is)...);
}

template <size_t N, typename F>
auto gen_tuple (const F & func) {
    return gen_tuple_impl(func, std::make_index_sequence<N>{} );
}




template<size_t Is = 0, typename... Tp>
void matrix_multiplication (std::vector<std::tuple<Tp...>> & A, std::vector<std::tuple<Tp...>> & B, std::vector<std::tuple<Tp...>> & R) {
    std::tuple<Tp...> column;
    std::vector<double> buf_column;
    buf_column.reserve(B.size());
    for (const auto & vector : B)
        buf_column.emplace_back(std::get<Is>(vector));
    column = std::move(gen_tuple<dim>([&](size_t j) { return buf_column[j]; }));
    for (int i = 0; i < B.size(); ++i)
        std::get<Is>(R[i]) = scalar_vector_multiplication(A[i], column);
    if constexpr (Is + 1 != sizeof...(Tp))
        matrix_multiplication<Is + 1>(A, B, R);
}


