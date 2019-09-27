// Generated by rstantools.  Do not edit by hand.

/*
    rstanlm is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    rstanlm is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with rstanlm.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef MODELS_HPP
#define MODELS_HPP
#define STAN__SERVICES__COMMAND_HPP
#include <rstan/rstaninc.hpp>
// Code generated by Stan version 2.19.1
#include <stan/model/model_header.hpp>
namespace model_sn_namespace {
using std::istream;
using std::string;
using std::stringstream;
using std::vector;
using stan::io::dump;
using stan::math::lgamma;
using stan::model::prob_grad;
using namespace stan::math;
static int current_statement_begin__;
stan::io::program_reader prog_reader__() {
    stan::io::program_reader reader;
    reader.add_event(0, 0, "start", "model_sn");
    reader.add_event(35, 33, "end", "model_sn");
    return reader;
}
#include <stan_meta_header.hpp>
class model_sn : public prob_grad {
private:
        int K;
        std::vector<double> yi;
        std::vector<double> si;
public:
    model_sn(stan::io::var_context& context__,
        std::ostream* pstream__ = 0)
        : prob_grad(0) {
        ctor_body(context__, 0, pstream__);
    }
    model_sn(stan::io::var_context& context__,
        unsigned int random_seed__,
        std::ostream* pstream__ = 0)
        : prob_grad(0) {
        ctor_body(context__, random_seed__, pstream__);
    }
    void ctor_body(stan::io::var_context& context__,
                   unsigned int random_seed__,
                   std::ostream* pstream__) {
        typedef double local_scalar_t__;
        boost::ecuyer1988 base_rng__ =
          stan::services::util::create_rng(random_seed__, 0);
        (void) base_rng__;  // suppress unused var warning
        current_statement_begin__ = -1;
        static const char* function__ = "model_sn_namespace::model_sn";
        (void) function__;  // dummy to suppress unused var warning
        size_t pos__;
        (void) pos__;  // dummy to suppress unused var warning
        std::vector<int> vals_i__;
        std::vector<double> vals_r__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        try {
            // initialize data block variables from context__
            current_statement_begin__ = 3;
            context__.validate_dims("data initialization", "K", "int", context__.to_vec());
            K = int(0);
            vals_i__ = context__.vals_i("K");
            pos__ = 0;
            K = vals_i__[pos__++];
            check_greater_or_equal(function__, "K", K, 0);
            current_statement_begin__ = 4;
            validate_non_negative_index("yi", "K", K);
            context__.validate_dims("data initialization", "yi", "double", context__.to_vec(K));
            yi = std::vector<double>(K, double(0));
            vals_r__ = context__.vals_r("yi");
            pos__ = 0;
            size_t yi_k_0_max__ = K;
            for (size_t k_0__ = 0; k_0__ < yi_k_0_max__; ++k_0__) {
                yi[k_0__] = vals_r__[pos__++];
            }
            current_statement_begin__ = 5;
            validate_non_negative_index("si", "K", K);
            context__.validate_dims("data initialization", "si", "double", context__.to_vec(K));
            si = std::vector<double>(K, double(0));
            vals_r__ = context__.vals_r("si");
            pos__ = 0;
            size_t si_k_0_max__ = K;
            for (size_t k_0__ = 0; k_0__ < si_k_0_max__; ++k_0__) {
                si[k_0__] = vals_r__[pos__++];
            }
            size_t si_i_0_max__ = K;
            for (size_t i_0__ = 0; i_0__ < si_i_0_max__; ++i_0__) {
                check_greater_or_equal(function__, "si[i_0__]", si[i_0__], 0);
            }
            // initialize transformed data variables
            // execute transformed data statements
            // validate transformed data
            // validate, set parameter ranges
            num_params_r__ = 0U;
            param_ranges_i__.clear();
            current_statement_begin__ = 8;
            num_params_r__ += 1;
            current_statement_begin__ = 9;
            num_params_r__ += 1;
            current_statement_begin__ = 10;
            num_params_r__ += 1;
            current_statement_begin__ = 11;
            validate_non_negative_index("theta", "K", K);
            num_params_r__ += (1 * K);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    ~model_sn() { }
    void transform_inits(const stan::io::var_context& context__,
                         std::vector<int>& params_i__,
                         std::vector<double>& params_r__,
                         std::ostream* pstream__) const {
        typedef double local_scalar_t__;
        stan::io::writer<double> writer__(params_r__, params_i__);
        size_t pos__;
        (void) pos__; // dummy call to supress warning
        std::vector<double> vals_r__;
        std::vector<int> vals_i__;
        current_statement_begin__ = 8;
        if (!(context__.contains_r("xi")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable xi missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("xi");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "xi", "double", context__.to_vec());
        double xi(0);
        xi = vals_r__[pos__++];
        try {
            writer__.scalar_unconstrain(xi);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable xi: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 9;
        if (!(context__.contains_r("omega")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable omega missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("omega");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "omega", "double", context__.to_vec());
        double omega(0);
        omega = vals_r__[pos__++];
        try {
            writer__.scalar_lb_unconstrain(0, omega);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable omega: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 10;
        if (!(context__.contains_r("alpha")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable alpha missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("alpha");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "alpha", "double", context__.to_vec());
        double alpha(0);
        alpha = vals_r__[pos__++];
        try {
            writer__.scalar_unconstrain(alpha);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable alpha: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 11;
        if (!(context__.contains_r("theta")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable theta missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("theta");
        pos__ = 0U;
        validate_non_negative_index("theta", "K", K);
        context__.validate_dims("parameter initialization", "theta", "double", context__.to_vec(K));
        std::vector<double> theta(K, double(0));
        size_t theta_k_0_max__ = K;
        for (size_t k_0__ = 0; k_0__ < theta_k_0_max__; ++k_0__) {
            theta[k_0__] = vals_r__[pos__++];
        }
        size_t theta_i_0_max__ = K;
        for (size_t i_0__ = 0; i_0__ < theta_i_0_max__; ++i_0__) {
            try {
                writer__.scalar_unconstrain(theta[i_0__]);
            } catch (const std::exception& e) {
                stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable theta: ") + e.what()), current_statement_begin__, prog_reader__());
            }
        }
        params_r__ = writer__.data_r();
        params_i__ = writer__.data_i();
    }
    void transform_inits(const stan::io::var_context& context,
                         Eigen::Matrix<double, Eigen::Dynamic, 1>& params_r,
                         std::ostream* pstream__) const {
      std::vector<double> params_r_vec;
      std::vector<int> params_i_vec;
      transform_inits(context, params_i_vec, params_r_vec, pstream__);
      params_r.resize(params_r_vec.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r(i) = params_r_vec[i];
    }
    template <bool propto__, bool jacobian__, typename T__>
    T__ log_prob(std::vector<T__>& params_r__,
                 std::vector<int>& params_i__,
                 std::ostream* pstream__ = 0) const {
        typedef T__ local_scalar_t__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // dummy to suppress unused var warning
        T__ lp__(0.0);
        stan::math::accumulator<T__> lp_accum__;
        try {
            stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
            // model parameters
            current_statement_begin__ = 8;
            local_scalar_t__ xi;
            (void) xi;  // dummy to suppress unused var warning
            if (jacobian__)
                xi = in__.scalar_constrain(lp__);
            else
                xi = in__.scalar_constrain();
            current_statement_begin__ = 9;
            local_scalar_t__ omega;
            (void) omega;  // dummy to suppress unused var warning
            if (jacobian__)
                omega = in__.scalar_lb_constrain(0, lp__);
            else
                omega = in__.scalar_lb_constrain(0);
            current_statement_begin__ = 10;
            local_scalar_t__ alpha;
            (void) alpha;  // dummy to suppress unused var warning
            if (jacobian__)
                alpha = in__.scalar_constrain(lp__);
            else
                alpha = in__.scalar_constrain();
            current_statement_begin__ = 11;
            std::vector<local_scalar_t__> theta;
            size_t theta_d_0_max__ = K;
            theta.reserve(theta_d_0_max__);
            for (size_t d_0__ = 0; d_0__ < theta_d_0_max__; ++d_0__) {
                if (jacobian__)
                    theta.push_back(in__.scalar_constrain(lp__));
                else
                    theta.push_back(in__.scalar_constrain());
            }
            // transformed parameters
            current_statement_begin__ = 14;
            local_scalar_t__ mu;
            (void) mu;  // dummy to suppress unused var warning
            stan::math::initialize(mu, DUMMY_VAR__);
            stan::math::fill(mu, DUMMY_VAR__);
            current_statement_begin__ = 15;
            local_scalar_t__ V;
            (void) V;  // dummy to suppress unused var warning
            stan::math::initialize(V, DUMMY_VAR__);
            stan::math::fill(V, DUMMY_VAR__);
            // transformed parameters block statements
            current_statement_begin__ = 16;
            stan::math::assign(mu, (xi + (((omega * stan::math::sqrt((2 / stan::math::pi()))) * alpha) / stan::math::sqrt((1 + pow(alpha, 2))))));
            current_statement_begin__ = 17;
            stan::math::assign(V, (pow(omega, 2) * (1 - pow(((stan::math::sqrt((2 / stan::math::pi())) * alpha) / stan::math::sqrt((1 + pow(alpha, 2)))), 2))));
            // validate transformed parameters
            const char* function__ = "validate transformed params";
            (void) function__;  // dummy to suppress unused var warning
            current_statement_begin__ = 14;
            if (stan::math::is_uninitialized(mu)) {
                std::stringstream msg__;
                msg__ << "Undefined transformed parameter: mu";
                stan::lang::rethrow_located(std::runtime_error(std::string("Error initializing variable mu: ") + msg__.str()), current_statement_begin__, prog_reader__());
            }
            current_statement_begin__ = 15;
            if (stan::math::is_uninitialized(V)) {
                std::stringstream msg__;
                msg__ << "Undefined transformed parameter: V";
                stan::lang::rethrow_located(std::runtime_error(std::string("Error initializing variable V: ") + msg__.str()), current_statement_begin__, prog_reader__());
            }
            check_greater_or_equal(function__, "V", V, 0);
            // model body
            current_statement_begin__ = 20;
            for (int i = 1; i <= K; ++i) {
                current_statement_begin__ = 21;
                lp_accum__.add(normal_log<propto__>(get_base1(yi, i, "yi", 1), get_base1(theta, i, "theta", 1), get_base1(si, i, "si", 1)));
                current_statement_begin__ = 22;
                lp_accum__.add(skew_normal_log<propto__>(get_base1(theta, i, "theta", 1), xi, omega, alpha));
            }
            current_statement_begin__ = 24;
            lp_accum__.add(normal_log<propto__>(xi, 0, 100));
            current_statement_begin__ = 25;
            lp_accum__.add(uniform_log<propto__>(omega, 0, 20));
            current_statement_begin__ = 26;
            lp_accum__.add(normal_log<propto__>(alpha, 0, 5));
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
        lp_accum__.add(lp__);
        return lp_accum__.sum();
    } // log_prob()
    template <bool propto, bool jacobian, typename T_>
    T_ log_prob(Eigen::Matrix<T_,Eigen::Dynamic,1>& params_r,
               std::ostream* pstream = 0) const {
      std::vector<T_> vec_params_r;
      vec_params_r.reserve(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        vec_params_r.push_back(params_r(i));
      std::vector<int> vec_params_i;
      return log_prob<propto,jacobian,T_>(vec_params_r, vec_params_i, pstream);
    }
    void get_param_names(std::vector<std::string>& names__) const {
        names__.resize(0);
        names__.push_back("xi");
        names__.push_back("omega");
        names__.push_back("alpha");
        names__.push_back("theta");
        names__.push_back("mu");
        names__.push_back("V");
        names__.push_back("log_lik");
        names__.push_back("theta_new");
    }
    void get_dims(std::vector<std::vector<size_t> >& dimss__) const {
        dimss__.resize(0);
        std::vector<size_t> dims__;
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(K);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(K);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
    }
    template <typename RNG>
    void write_array(RNG& base_rng__,
                     std::vector<double>& params_r__,
                     std::vector<int>& params_i__,
                     std::vector<double>& vars__,
                     bool include_tparams__ = true,
                     bool include_gqs__ = true,
                     std::ostream* pstream__ = 0) const {
        typedef double local_scalar_t__;
        vars__.resize(0);
        stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
        static const char* function__ = "model_sn_namespace::write_array";
        (void) function__;  // dummy to suppress unused var warning
        // read-transform, write parameters
        double xi = in__.scalar_constrain();
        vars__.push_back(xi);
        double omega = in__.scalar_lb_constrain(0);
        vars__.push_back(omega);
        double alpha = in__.scalar_constrain();
        vars__.push_back(alpha);
        std::vector<double> theta;
        size_t theta_d_0_max__ = K;
        theta.reserve(theta_d_0_max__);
        for (size_t d_0__ = 0; d_0__ < theta_d_0_max__; ++d_0__) {
            theta.push_back(in__.scalar_constrain());
        }
        size_t theta_k_0_max__ = K;
        for (size_t k_0__ = 0; k_0__ < theta_k_0_max__; ++k_0__) {
            vars__.push_back(theta[k_0__]);
        }
        double lp__ = 0.0;
        (void) lp__;  // dummy to suppress unused var warning
        stan::math::accumulator<double> lp_accum__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        if (!include_tparams__ && !include_gqs__) return;
        try {
            // declare and define transformed parameters
            current_statement_begin__ = 14;
            double mu;
            (void) mu;  // dummy to suppress unused var warning
            stan::math::initialize(mu, DUMMY_VAR__);
            stan::math::fill(mu, DUMMY_VAR__);
            current_statement_begin__ = 15;
            double V;
            (void) V;  // dummy to suppress unused var warning
            stan::math::initialize(V, DUMMY_VAR__);
            stan::math::fill(V, DUMMY_VAR__);
            // do transformed parameters statements
            current_statement_begin__ = 16;
            stan::math::assign(mu, (xi + (((omega * stan::math::sqrt((2 / stan::math::pi()))) * alpha) / stan::math::sqrt((1 + pow(alpha, 2))))));
            current_statement_begin__ = 17;
            stan::math::assign(V, (pow(omega, 2) * (1 - pow(((stan::math::sqrt((2 / stan::math::pi())) * alpha) / stan::math::sqrt((1 + pow(alpha, 2)))), 2))));
            if (!include_gqs__ && !include_tparams__) return;
            // validate transformed parameters
            const char* function__ = "validate transformed params";
            (void) function__;  // dummy to suppress unused var warning
            current_statement_begin__ = 15;
            check_greater_or_equal(function__, "V", V, 0);
            // write transformed parameters
            if (include_tparams__) {
                vars__.push_back(mu);
                vars__.push_back(V);
            }
            if (!include_gqs__) return;
            // declare and define generated quantities
            current_statement_begin__ = 29;
            validate_non_negative_index("log_lik", "K", K);
            Eigen::Matrix<double, Eigen::Dynamic, 1> log_lik(K);
            stan::math::initialize(log_lik, DUMMY_VAR__);
            stan::math::fill(log_lik, DUMMY_VAR__);
            current_statement_begin__ = 30;
            double theta_new;
            (void) theta_new;  // dummy to suppress unused var warning
            stan::math::initialize(theta_new, DUMMY_VAR__);
            stan::math::fill(theta_new, DUMMY_VAR__);
            // generated quantities statements
            current_statement_begin__ = 31;
            for (int i = 1; i <= K; ++i) {
                current_statement_begin__ = 31;
                stan::model::assign(log_lik, 
                            stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list()), 
                            normal_log(get_base1(yi, i, "yi", 1), get_base1(theta, i, "theta", 1), get_base1(si, i, "si", 1)), 
                            "assigning variable log_lik");
            }
            current_statement_begin__ = 32;
            stan::math::assign(theta_new, skew_normal_rng(xi, omega, alpha, base_rng__));
            // validate, write generated quantities
            current_statement_begin__ = 29;
            size_t log_lik_j_1_max__ = K;
            for (size_t j_1__ = 0; j_1__ < log_lik_j_1_max__; ++j_1__) {
                vars__.push_back(log_lik(j_1__));
            }
            current_statement_begin__ = 30;
            vars__.push_back(theta_new);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    template <typename RNG>
    void write_array(RNG& base_rng,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& vars,
                     bool include_tparams = true,
                     bool include_gqs = true,
                     std::ostream* pstream = 0) const {
      std::vector<double> params_r_vec(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r_vec[i] = params_r(i);
      std::vector<double> vars_vec;
      std::vector<int> params_i_vec;
      write_array(base_rng, params_r_vec, params_i_vec, vars_vec, include_tparams, include_gqs, pstream);
      vars.resize(vars_vec.size());
      for (int i = 0; i < vars.size(); ++i)
        vars(i) = vars_vec[i];
    }
    static std::string model_name() {
        return "model_sn";
    }
    void constrained_param_names(std::vector<std::string>& param_names__,
                                 bool include_tparams__ = true,
                                 bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        param_name_stream__.str(std::string());
        param_name_stream__ << "xi";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "omega";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "alpha";
        param_names__.push_back(param_name_stream__.str());
        size_t theta_k_0_max__ = K;
        for (size_t k_0__ = 0; k_0__ < theta_k_0_max__; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "theta" << '.' << k_0__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "mu";
            param_names__.push_back(param_name_stream__.str());
            param_name_stream__.str(std::string());
            param_name_stream__ << "V";
            param_names__.push_back(param_name_stream__.str());
        }
        if (!include_gqs__) return;
        size_t log_lik_j_1_max__ = K;
        for (size_t j_1__ = 0; j_1__ < log_lik_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "log_lik" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        param_name_stream__.str(std::string());
        param_name_stream__ << "theta_new";
        param_names__.push_back(param_name_stream__.str());
    }
    void unconstrained_param_names(std::vector<std::string>& param_names__,
                                   bool include_tparams__ = true,
                                   bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        param_name_stream__.str(std::string());
        param_name_stream__ << "xi";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "omega";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "alpha";
        param_names__.push_back(param_name_stream__.str());
        size_t theta_k_0_max__ = K;
        for (size_t k_0__ = 0; k_0__ < theta_k_0_max__; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "theta" << '.' << k_0__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "mu";
            param_names__.push_back(param_name_stream__.str());
            param_name_stream__.str(std::string());
            param_name_stream__ << "V";
            param_names__.push_back(param_name_stream__.str());
        }
        if (!include_gqs__) return;
        size_t log_lik_j_1_max__ = K;
        for (size_t j_1__ = 0; j_1__ < log_lik_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "log_lik" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        param_name_stream__.str(std::string());
        param_name_stream__ << "theta_new";
        param_names__.push_back(param_name_stream__.str());
    }
}; // model
}  // namespace
typedef model_sn_namespace::model_sn stan_model;
#endif
