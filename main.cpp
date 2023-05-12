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

typedef std::vector<point> matrix;
typedef std::vector<matrix> cation_frame;
typedef std::vector<cation_frame> cation_frames;

typedef std::vector<point> frame;
typedef std::vector<frame> frames;


const size_t dim = std::tuple_size<point>{}; // Number of spatial dimensions.


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


template<typename... Tp>
void row_addition (std::vector<std::tuple<Tp...>> & matrix, double & lambda, const size_t & i, const size_t & j) {
    std::tuple<Tp...> R_i;
    vector_scalar_multiplication(matrix[i], lambda, R_i);
    vector_offset(matrix[j], R_i, matrix[j]);
}


template<typename... Tp>
void row_switching (std::vector<std::tuple<Tp...>> & matrix, const size_t & i, const size_t & j) {
    std::swap(matrix[i], matrix[j]);
    if (j%2 == 0)
        vector_scalar_multiplication(matrix[1], -1, matrix[1]);
}


bool is_equal(const double & x, const double & y) {
    return std::fabs(x - y) < std::numeric_limits<double>::epsilon();
}


template<size_t Is = 0, typename... Tp>
void forward (std::vector<std::tuple<Tp...>> & invertible_matrix, std::vector<std::tuple<Tp...>> & identity_matrix) {
    for (int i = Is+1; i < invertible_matrix.size(); ++i) {
        int k = 0;
        while (is_equal(std::get<Is>(invertible_matrix[Is]), 0) && k < invertible_matrix.size()) {
            row_switching(invertible_matrix,Is, k);
            row_switching(identity_matrix, Is, k);
            ++k;
        }
        double lambda = -std::get<Is>(invertible_matrix[i]) / std::get<Is>(invertible_matrix[Is]);
        row_addition(invertible_matrix, lambda, Is, i);
        row_addition(identity_matrix, lambda, Is, i);
    }
    if constexpr (Is + 1 != sizeof...(Tp))
        forward<Is + 1>(invertible_matrix, identity_matrix);
}


template<size_t Is = 0, typename... Tp, size_t N = sizeof...(Tp)-1>
void back (std::vector<std::tuple<Tp...>> & invertible_matrix, std::vector<std::tuple<Tp...>> & identity_matrix) {
    for (int i = N-Is; i > 0; --i) {
        int k = 0;
        while (is_equal(std::get<N-Is>(invertible_matrix[N-Is]), 0) && k < invertible_matrix.size()) {
            row_switching(invertible_matrix, N-Is, k);
            row_switching(identity_matrix, N-Is, k);
            ++k;
        }
        double lambda = -std::get<N-Is>(invertible_matrix[i-1]) / std::get<N-Is>(invertible_matrix[N-Is]);
        row_addition(invertible_matrix, lambda, N-Is, i-1);
        row_addition(identity_matrix, lambda, N-Is, i-1);
    }
    if constexpr (N-Is > 0)
        back<Is + 1>(invertible_matrix, identity_matrix);
}


template<size_t Is = 0, typename... Tp>
void final_step (std::vector<std::tuple<Tp...>> & invertible_matrix, std::vector<std::tuple<Tp...>> & identity_matrix) {
    vector_scalar_multiplication(identity_matrix[Is], 1.0 / std::get<Is>(invertible_matrix[Is]), identity_matrix[Is]);
    for (int i = 0; i < identity_matrix.size(); ++i)
        if(is_equal(std::get<Is>(identity_matrix[i]), -0))
            std::get<Is>(identity_matrix[i]) = 0; // Just for output due to DRY.
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


template<typename... Tp>
auto Gaussian_elimination (std::vector<std::tuple<Tp...>> invertible_matrix) { // !&?
    std::vector<std::tuple<Tp...>> E (invertible_matrix.size());
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
void matrix_multiplication (const std::vector<std::tuple<Tp...>> & A, const std::vector<std::tuple<Tp...>> & B, std::vector<std::tuple<Tp...>> & R) {
    std::tuple<Tp...> column;
    std::vector<double> buf_column;
    buf_column.reserve(B.size());
    for (const auto & vector : B)
        buf_column.emplace_back(std::get<Is>(vector));
    column = std::move(gen_tuple<dim>([&](size_t j) { return buf_column[j]; }));
    for (int i = 0; i < B.size(); ++i)
        std::get<Is>(R[i]) = std::move(scalar_vector_multiplication(A[i], column));
    if constexpr (Is + 1 != sizeof...(Tp))
        matrix_multiplication<Is + 1>(A, B, R);
}


template<size_t Is = 0, typename... Tp>
void matrix_transposition (const std::vector<std::tuple<Tp...>> & A, std::vector<std::tuple<Tp...>> & R) {
    std::vector<double> buf_column;
    buf_column.reserve(A.size());
    for (const auto & vector : A)
        buf_column.emplace_back(std::get<Is>(vector));
    R[Is] = std::move(gen_tuple<dim>([&](size_t j) { return buf_column[j]; }));
    if constexpr (Is + 1 != sizeof...(Tp)) // !
        matrix_transposition<Is + 1>(A,  R);
}


template<size_t Is = 0, typename... Tp>
void get_column (const std::vector<std::tuple<Tp...>> & A, const int & column_number, std::tuple<Tp...> & column) {
    if (Is == column_number) {
        std::vector<double> buf_column;
        buf_column.reserve(A.size());
        for (const auto & vector: A)
            buf_column.emplace_back(std::get<Is>(vector));
        column = std::move(gen_tuple<dim>([&](size_t j) { return buf_column[j]; }));
    }
    if constexpr (Is + 1 != sizeof...(Tp))
        get_column<Is + 1>(A,  column_number, column);
}


void direction (double & cos_a, double & cos_b, double & cos_g, double & sin_a, double & sin_b, double & sin_g,
                point & direction_vector, const std::vector<point> & basis_set) {

    point xOy_projection = std::make_tuple(std::get<0>(direction_vector), std::get<1>(direction_vector), 0);

    cos_a = cos_angle_between_vectors(xOy_projection, basis_set[0]);
    sin_a = sin_angle_between_vectors(xOy_projection, basis_set[0]);

    cos_b = cos_angle_between_vectors(xOy_projection, basis_set[1]);
    sin_b = sin_angle_between_vectors(xOy_projection, basis_set[1]);

    cos_g = cos_angle_between_vectors(direction_vector, basis_set[2]);
    sin_g = sin_angle_between_vectors(direction_vector, basis_set[2]);
}


std::vector<point> frame_of_reference_rotation (point & direction_vector, const std::vector<point> & basis_set) {
    std::vector<point> new_basis_set (basis_set.size());

    double cos_a, cos_b, cos_g, sin_a, sin_b, sin_g;

    direction (cos_a, cos_b, cos_g, sin_a, sin_b, sin_g, direction_vector, basis_set);

    std::vector<point> M; // Rotation matrix
    M.emplace_back(cos_b*cos_g,                     -sin_a*cos_b,                      sin_b);
    M.emplace_back(sin_a*sin_b*cos_g + sin_g*cos_a, -sin_a*sin_b*sin_g + cos_a*cos_g, -sin_a*cos_b);
    M.emplace_back(sin_a*sin_g - sin_b*cos_a*cos_g,  sin_a*cos_g + sin_b*sin_g*cos_a,  cos_a*cos_b);

    auto Revers_rotation_matrix = std::move(Gaussian_elimination(M));
    matrix_multiplication(Revers_rotation_matrix, basis_set, new_basis_set);

    return new_basis_set;
}


point old_vector_in_new_basis_set (point & vector, std::vector<point> & basis_set) {
    point result;
    std::vector<point> completed (basis_set.size());
    completed[0] = vector; // May be matrix_insert_template?
    matrix_transposition(completed, completed);
    matrix_multiplication(basis_set, completed, completed);
    get_column(completed, 0, result);
    return result;
}


template<typename T, size_t... Is>
bool outside_the_box_impl (T const& t, std::index_sequence<Is...>, const double & left_border, const double & right_border) {
    return ((std::get<Is>(t) > right_border || std::get<Is>(t) < left_border) | ...);
}

template <class Tuple>
bool outside_the_box (const Tuple& t, const double & left_border, const double & right_border) {
    constexpr auto size = std::tuple_size<Tuple>{};
    return outside_the_box_impl(t, std::make_index_sequence<size>{}, left_border, right_border);
}


void outsiders_excluding (frame & atoms, const double & left_border, const double & right_border) {
    for (int j = 0; j < atoms.size(); ++j)
        if (outside_the_box(atoms[j], left_border, right_border))
            atoms.erase(atoms.begin()+j-1);
}


void data_clearing (frames & oxygens, frames & free_protons, const double & left_border, const double & right_border) {
    for (long i = 0; i < oxygens.size(); ++i) {
        outsiders_excluding(oxygens[i], left_border, right_border);
        outsiders_excluding(free_protons[i], left_border, right_border);
    }
}


bool cation_sufficient_condition (point & O1, point & O2, point & H) {
    double dist = distance(O1, O2);
    double dist21st = distance(H, O1);
    double dist22nd = distance(H, O2);
    return ((dist > dist21st && dist > dist22nd));
}


std::vector<point> oxygen_neighbour_list (int & oxygen_index, frame & oxygens, const double & length) {
    std::vector<point> list;
    for (int i = 0; i < oxygens.size() ; ++i)
        if (distance(oxygens[i], oxygens[oxygen_index]) <= length)
            list.emplace_back(oxygens[i]);
    return list;
}


template<size_t Is = 0, typename... Tp>
void diag_multiplication (std::vector<std::tuple<Tp...>> & M, double & result) {
    result *= std::get<Is>(M[Is]);
    if constexpr (Is + 1 != sizeof...(Tp))
        diag_multiplication<Is + 1>(M, result);
}


template<typename... Tp>
double det (std::vector<std::tuple<Tp...>>  M) {
    double result = 1;
    std::vector<std::tuple<Tp...>> foo (M.size());
    forward(M, foo);
    diag_multiplication(M, result);
    return result;
}





//cation_frame is_cation (frame & oxygens, frame & free_protons, const double & deviation,
//                        const double & left_border, const double & right_border) {
//
//
//
//
//    std::vector<int> oxygens_excludes, protons_excludes;
//
//    for (int i = 0; i < oxygens.size(); ++i)
//        if (outside_the_box(oxygens[i], left_border, right_border))
//            oxygens_excludes.emplace_back(i);
//
//
//    for (int i = 0; i < oxygens.size(); ++i) {
//
//    }
//
//
//
//
//    // Если катион, то необходимо удалить из фреймов соответсвующие частицы. Тогда они должны пода
//    // Можно удалять сразу окружение по критерию.
//
//
//    for (int i = 0; i < oxygens.size(); ++i) { // add ! any_of excludes
//
//        if (std::any_of(oxygens_excludes.begin(), oxygens_excludes.end(), [&](int k) { return k == i; }))
//            continue;
//
//        point ref, O_2_in_O1_ref; // здесь критерий нужен
//        vector_scalar_multiplication(O_1, -1.0, ref); // !
//        vector_offset(O_2, ref, O_2_in_O1_ref);
//
//        const std::vector<point> default_basis_set = {std::make_tuple(1, 0, 0),
//                                                      std::make_tuple(0, 1, 0),
//                                                      std::make_tuple(0, 0, 1)};
//        // You lost offset
//        std::vector<point> new_basis_set = std::move(frame_of_reference_rotation(O_2_in_O1_ref, default_basis_set));
//
//
//        for (int i = 0; i < free_protons.size(); ++i) {
//            vector_offset(, ref, O_2_in_O1_ref);
//        }
//    }
//    return ans;
//}
//
//
//bool pair_of_interest (point & O_1, point & O_2, frame & free_proton, const double & deviation, const std::string & type) {
//    if (type == "Zundel") {
//        point geometric_center = geometric_center_between_O(O_1,O_2);
//    }
//
//}


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


