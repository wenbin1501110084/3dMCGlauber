// Copyright @ Chun Shen 2018

#include "Glauber.h"
#include "data_structs.h"
#include "PhysConsts.h"

#include <string>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <iostream>
#include <algorithm>
#include <memory>

using std::cout;
using std::endl;
using std::shared_ptr;

// put in use of baryon used flags

namespace MCGlb {

Glauber::Glauber(const MCGlb::Parameters &param_in,
                 shared_ptr<RandomUtil::Random> ran_gen) :
    parameter_list(param_in) {
    sample_valence_quark = false;
    if (parameter_list.get_use_quarks() > 0) {
        sample_valence_quark = true;
    }
    projectile = std::unique_ptr<Nucleus>(
            new Nucleus(parameter_list.get_projectle_nucleus_name(), ran_gen,
                        sample_valence_quark));
    target = std::unique_ptr<Nucleus>(
            new Nucleus(parameter_list.get_target_nucleus_name(), ran_gen,
                        sample_valence_quark));
    if (sample_valence_quark) {
        projectile->set_valence_quark_Q2(parameter_list.get_quarks_Q2());
        target->set_valence_quark_Q2(parameter_list.get_quarks_Q2());
    }
    ran_gen_ptr = ran_gen;
    string_production_mode = parameter_list.get_QCD_string_production_mode();
}

void Glauber::make_nuclei() {
    projectile->generate_nucleus_3d_configuration();
    target->generate_nucleus_3d_configuration();
    projectile->accelerate_nucleus(parameter_list.get_roots(), 1);
    target->accelerate_nucleus(parameter_list.get_roots(), -1);

    // sample impact parameters
    auto b_max = parameter_list.get_b_max();
    auto b_min = parameter_list.get_b_min();
    impact_b = sqrt(b_min*b_min +
            (b_max*b_max - b_min*b_min)*ran_gen_ptr.lock()->rand_uniform());
    SpatialVec proj_shift = {0., impact_b/2., 0.,
                             -projectile->get_z_max() - 1e-15};
    projectile->shift_nucleus(proj_shift);
    SpatialVec targ_shift = {0., -impact_b/2., 0.,
                             -target->get_z_min() + 1e-15};
    target->shift_nucleus(targ_shift);
    //projectile->output_nucleon_positions("projectile.dat");
    //target->output_nucleon_positions("target.dat");
}


int Glauber::make_collision_schedule() {
    collision_schedule.clear();
    auto d2 = (compute_NN_inelastic_cross_section(parameter_list.get_roots())
               /(M_PI*10.));  // in fm^2 
    auto proj_nucleon_list = projectile->get_nucleon_list();
    auto targ_nucleon_list = target->get_nucleon_list();
    for (auto &iproj: (*proj_nucleon_list)) {
        auto proj_x = iproj->get_x();
        for (auto &itarg: (*targ_nucleon_list)) {
            auto targ_x = itarg->get_x();
            auto dij = (  (targ_x[1] - proj_x[1])*(targ_x[1] - proj_x[1])
                        + (targ_x[2] - proj_x[2])*(targ_x[2] - proj_x[2]));
            if (hit(dij, d2)) {
                create_a_collision_event(iproj, itarg);
                iproj->set_wounded(true);
                itarg->set_wounded(true);
                iproj->add_collide_nucleon(itarg);
                itarg->add_collide_nucleon(iproj);
            }
        }
    }
    return(collision_schedule.size());
}

int Glauber::get_Npart() const {
    int Npart = (projectile->get_number_of_wounded_nucleons()
                  + target->get_number_of_wounded_nucleons());
    return(Npart);
}

bool Glauber::hit(real d2, real d2_in) const {
    real G = 0.92;  // from Glassando
    return(ran_gen_ptr.lock()->rand_uniform() < G*exp(-G*d2/d2_in));
}

void Glauber::create_a_collision_event(weak_ptr<Nucleon> proj,
                                       weak_ptr<Nucleon> targ) {
    real t_coll, z_coll;
    auto x1 = proj.lock()->get_x();
    auto p1 = proj.lock()->get_p();
    auto x2 = targ.lock()->get_x();
    auto p2 = targ.lock()->get_p();
    auto v1 = p1[3]/p1[0];
    auto v2 = p2[3]/p2[0];
    auto collided = get_collision_point(x1[0], x1[3], v1, x2[3], v2,
                                        t_coll, z_coll);
    if (collided) {
        SpatialVec x_coll = {t_coll, (x1[1] + x2[1])/2.,
                             (x1[2] + x2[2])/2., z_coll};
        CollisionEvent event_ptr(x_coll, proj, targ);
        if (proj.lock()->is_connected_with(targ)) {
            event_ptr.set_produced_a_string(true);
        }
        collision_schedule.insert(event_ptr);
    }
}

bool Glauber::get_collision_point(real t, real z1, real v1, real z2, real v2,
                                  real &t_coll, real &z_coll) const {
    bool collided = false;
    real delta_t = (z2 - z1)/(v1 - v2);
    if (delta_t > 0.) {
        collided = true;
        t_coll = t  + delta_t;
        z_coll = z1 + v1*delta_t;
    } else {
        t_coll = -1.;
        z_coll = -1.;
    }
    return(collided);
}

real Glauber::compute_NN_inelastic_cross_section(real ecm) const {
    real s = ecm*ecm;
    real sigma_NN_total = 44.4 - 2.9*log(s) + 0.33*log(s)*log(s);
    real sigma_NN_inel  = (sigma_NN_total
                           - (11.4 - 1.52*log(s) + 0.13*log(s)*log(s)));
    return(sigma_NN_inel);
}


bool Glauber::decide_produce_string(CollisionEvent &event_ptr) const {
    bool flag_form_a_string = false;
    auto proj = event_ptr.get_proj_nucleon_ptr();
    auto targ = event_ptr.get_targ_nucleon_ptr();
    if (   proj.lock()->get_number_of_connections() == 0
        || targ.lock()->get_number_of_connections() == 0) {
      //bps: string gets flag for whether it has left or right or both baryon numbers.
      
      //int n_connects = (  proj->get_number_of_connections()
        //                  + targ->get_number_of_connections());
        //real R_A = pow(projectile->get_nucleus_A(), 1./3.);
        //real R_B = pow(target->get_nucleus_A(), 1./3.);
        //real asymmetry_factor = (2./std::max(R_A, R_B))*R_A*R_B/(R_A + R_B);
        //real cost_function = exp(- asymmetry_factor
        //                           *static_cast<real>(n_connects));
        //if (ran_gen_ptr.lock()->rand_uniform() < cost_function)
        //    flag_form_a_string = true;
        flag_form_a_string = true;
    }
    return(flag_form_a_string);
}


int Glauber::decide_QCD_strings_production() {
    std::vector<CollisionEvent> collision_list;
    for (auto &it: collision_schedule)
        collision_list.push_back(it);  // collision list is time ordered

    const auto QCD_string_production_mode =
                            parameter_list.get_QCD_string_production_mode();
    if (QCD_string_production_mode == 1) {
        // randomly ordered strings
        std::random_shuffle(collision_list.begin(), collision_list.end());
    } else if (QCD_string_production_mode == 2) {
        // anti-time ordered strings
        std::reverse(collision_list.begin(), collision_list.end());
    }

    int number_of_strings = 0;
    bool finished = false;
    while (!finished) {
        finished = true;
        for (auto &ievent: collision_list) {
            auto flag_form_a_string = decide_produce_string(ievent);
            auto proj = ievent.get_proj_nucleon_ptr();
            auto targ = ievent.get_targ_nucleon_ptr();
            if (flag_form_a_string) {
                number_of_strings++;
                proj.lock()->add_connected_nucleon(targ);
                targ.lock()->add_connected_nucleon(proj);
                ievent.set_produced_a_string(true);
            } else {
                if (   proj.lock()->get_number_of_connections() == 0
                    || targ.lock()->get_number_of_connections() == 0)
                        finished = false;
            }
        }
    }
    return(number_of_strings);
}


//! This function propagate individual nucleon inside the nucleus by dt
void Glauber::propagate_nuclei(real dt) {
    for (auto &n_i: (*projectile->get_nucleon_list()))
        propagate_nucleon(n_i, dt);
    for (auto &n_i: (*target->get_nucleon_list()))
        propagate_nucleon(n_i, dt);
}


void Glauber::propagate_nucleon(shared_ptr<Nucleon> n_i, real dt) {
    auto x_vec = n_i->get_x();
    auto p_vec = n_i->get_p();
    for (int i = 0; i < 4; i++) {
        real v_i = p_vec[i]/p_vec[0];
        x_vec[i] += dt*v_i;
    }
    n_i->set_x(x_vec);
}


void Glauber::update_momentum(shared_ptr<Nucleon> n_i, real y_shift) {
    auto pvec = n_i->get_p();
    auto y_i = atanh(pvec[3]/pvec[0]);
    auto y_f = y_i + y_shift;
    pvec[0] = PhysConsts::MProton*cosh(y_f);
    pvec[3] = PhysConsts::MProton*sinh(y_f);
    n_i->set_p(pvec);
}

// add flag if baryon used to proj and target
// - check it - only put in baryon if not used yet, then set it to used
int Glauber::perform_string_production() {
    if (sample_valence_quark) {
        projectile->sample_valence_quarks_inside_nucleons(
                                    parameter_list.get_roots(), 1);
        target->sample_valence_quarks_inside_nucleons(
                                    parameter_list.get_roots(), -1);
    }
    QCD_string_list.clear();
    const auto string_evolution_mode = (
                    parameter_list.get_QCD_string_evolution_mode());
    const auto baryon_junctions = parameter_list.get_baryon_junctions();
    bool has_baryon_left;
    bool has_baryon_right;
    real y_baryon_left;
    real y_baryon_right;
    
    // sqrt(parameter_list.get_roots()); // ~s^{-1/4} 
    real lambdaB =  parameter_list.get_lambdaB();
    lambdaB = std::min(1., lambdaB);
    
    //cout << lambdaB <<endl;

    real t_current = 0.0;
    int number_of_collided_events = 0;
    while (collision_schedule.size() > 0) {
        auto first_event = (*collision_schedule.begin());
        if (!first_event.is_valid()) {
            collision_schedule.erase((*collision_schedule.begin()));
            continue;
        }
        number_of_collided_events++;
        real dt = first_event.get_collision_time() - t_current;
        assert(dt > 0.);
        t_current = first_event.get_collision_time();
        propagate_nuclei(dt);
        if (!first_event.is_produced_a_string()) {
            collision_schedule.erase((*collision_schedule.begin()));
            continue;
        }
        auto x_coll = first_event.get_collision_position();
        auto proj = first_event.get_proj_nucleon_ptr().lock();
        auto targ = first_event.get_targ_nucleon_ptr().lock();
        real y_in_lrf = std::abs(proj->get_rapidity()
                                 - targ->get_rapidity())/2.;
        std::weak_ptr<Quark> proj_q;
        std::weak_ptr<Quark> targ_q;
        if (sample_valence_quark) {
            proj_q = proj->get_a_valence_quark();
            targ_q = targ->get_a_valence_quark();
            y_in_lrf = std::abs(proj_q.lock()->get_rapidity()
                                - targ_q.lock()->get_rapidity())/2.;
        }
        real tau_form = 0.5;      // [fm]
        real m_over_sigma = 1.0;  // [fm]
        if (string_evolution_mode == 1) {
            // fixed rapidity loss
            tau_form = 0.5;
            m_over_sigma = 1.0;
        } else if (string_evolution_mode == 2) {
            // both tau_form and sigma fluctuate
            auto y_loss = sample_rapidity_loss_shell(y_in_lrf);
            tau_form = 0.5 + 2.*ran_gen_ptr.lock()->rand_uniform();
            m_over_sigma = tau_form/sqrt(2.*(cosh(y_loss) - 1.));
        } else if (string_evolution_mode == 3) {
            // only tau_form fluctuates
            auto y_loss = sample_rapidity_loss_shell(y_in_lrf);
            tau_form = m_over_sigma*sqrt(2.*(cosh(y_loss) - 1.));
        } else if (string_evolution_mode == 4) {
            // only m_over_sigma fluctuates
            auto y_loss = sample_rapidity_loss_shell(y_in_lrf);
            m_over_sigma = tau_form/sqrt(2.*(cosh(y_loss) - 1.));
        }
        // set variables in case of no baryon junction transport
        has_baryon_left = 0;
        has_baryon_right = 0;
        y_baryon_right = 0.;
        y_baryon_left = 0.;
        if (baryon_junctions) {
            if (!proj->baryon_was_used()) {
                has_baryon_right = 1;
                proj->set_baryon_used(1);
            }
            if (!targ->baryon_was_used()) {
                has_baryon_left = 1;
                targ->set_baryon_used(1);
            } 
        }
        if (!sample_valence_quark) {
            QCDString qcd_string(x_coll, tau_form, proj, targ, m_over_sigma,
                                 has_baryon_right, has_baryon_left);
            QCD_string_list.push_back(qcd_string);
        } else {
            QCDString qcd_string(x_coll, tau_form, proj, targ,
                                 proj_q.lock(), targ_q.lock(), m_over_sigma,
                                 has_baryon_right, has_baryon_left);
            QCD_string_list.push_back(qcd_string);
        }
        real y_shift = 0.001;
        update_momentum(proj, -y_shift);
        update_momentum(targ,  y_shift);
        update_collision_schedule(first_event);
        collision_schedule.erase((*collision_schedule.begin()));
    }
    for (auto &it: QCD_string_list) {
        it.evolve_QCD_string();
        if (!baryon_junctions) {
            // set baryon rapidities to string endpoint rapidities
            // if no junction transport is used
            it.set_final_baryon_rapidities(it.get_y_f_left(),
                                           it.get_y_f_right());
        } else {
            // sample HERE if baryon should be moved
            y_baryon_right = 0.;
            y_baryon_left = 0.;
            if (it.get_has_baryon_right()) {
                if (ran_gen_ptr.lock()->rand_uniform() < lambdaB) {
                    // y_baryon_right = sample_junction_rapidity_right( it->get_y_i_left(), it->get_y_i_right() );
                    y_baryon_right = sample_junction_rapidity_right(
                                    it.get_y_f_left(), it.get_y_f_right());
                } else {
                    // One should use the very initial rapidities of the colliding nucleons to determine the junction rapidity.
                    // this may, however, lead to junctions lying outside the rapidity range of the final string which may caus problems
                    // in hydro and could be seen as unphysical...
                    // Another thing to discuss is the energy dependence. The cross section for baryon junctions stopping has a different root-s 
                    // dependence than the usual inelastic cross section - so in principle it needs to be treated separately altogether... not sure how yet
                    y_baryon_right = it.get_y_f_right();
                }
            }
            if (it.get_has_baryon_left()) {
                if (ran_gen_ptr.lock()->rand_uniform() < lambdaB) {
                    //y_baryon_left = sample_junction_rapidity_left( it->get_y_i_left(), it->get_y_i_right() );
                    y_baryon_left = sample_junction_rapidity_left(
                                    it.get_y_f_left(), it.get_y_f_right());
                } else {
                    y_baryon_left = it.get_y_f_left();
                }
            }
            it.set_final_baryon_rapidities(y_baryon_left,y_baryon_right);
        }
    }
    return(number_of_collided_events);
}


void Glauber::update_collision_schedule(CollisionEvent &event_happened) {
    auto proj = event_happened.get_proj_nucleon_ptr();
    proj.lock()->increment_collided_times();
    for (auto &it: (*proj.lock()->get_collide_nucleon_list()))
        create_a_collision_event(proj, it);
    auto targ = event_happened.get_targ_nucleon_ptr();
    targ.lock()->increment_collided_times();
    for (auto &it: (*targ.lock()->get_collide_nucleon_list()))
        create_a_collision_event(it, targ);
}

void Glauber::output_QCD_strings(std::string filename) {
    std::ofstream output(filename.c_str());
    output << "# norm  m_over_sigma[fm]  tau_form[fm]  tau_0[fm]  eta_s_0  "
           << "x_perp[fm]  y_perp[fm]  "
           << "eta_s_left  eta_s_right  y_l  y_r  fraction_l  fraction_r "
           << "y_l_i  y_r_i "
           << "eta_s_baryon_left  eta_s_baryon_right  y_l_baryon  y_r_baryon  "
           << endl;
    const auto baryon_junctions = (
                                   parameter_list.get_baryon_junctions());
    
    real fraction_left = 0.;
    real fraction_right = 0.;

    for (auto &it: QCD_string_list) {
        auto x_prod = it.get_x_production();
        auto tau_0  = sqrt(x_prod[0]*x_prod[0] - x_prod[3]*x_prod[3]);
        auto etas_0 = 0.5*log((x_prod[0] + x_prod[3])/(x_prod[0] - x_prod[3]));
        
        if (!baryon_junctions) {
            fraction_left = 1./(static_cast<real>(
                        it.get_proj().lock()->get_number_of_connections()));
            fraction_right = 1./(static_cast<real>(
                        it.get_targ().lock()->get_number_of_connections()));
        } else {
            fraction_left  = it.get_has_baryon_left();
            fraction_right = it.get_has_baryon_right();
        }

        real output_array[] = {
            1.0, it.get_m_over_sigma(), it.get_tau_form(),
            tau_0, etas_0, x_prod[1], x_prod[2],
            it.get_eta_s_left(), it.get_eta_s_right(),
            it.get_y_f_left(), it.get_y_f_right(),
            fraction_left, fraction_right,
            it.get_y_i_left(), it.get_y_i_right(),
            it.get_eta_s_baryon_left(), it.get_eta_s_baryon_right(),
            it.get_y_f_baryon_left(), it.get_y_f_baryon_right()
        };

        output << std::scientific << std::setprecision(8);
        for (int i = 0; i < 19; i++) {
            output << std::setw(15) << output_array[i] << "  ";
        }
        output << endl;
    }
    output.close();
}


real Glauber::sample_rapidity_loss_shell(real y_init) const {
    real y_loss = 0.0;
    if (parameter_list.get_rapidity_loss_method() == 1) {
        y_loss = sample_rapidity_loss_from_the_LEXUS_model(y_init);
    } else if (parameter_list.get_rapidity_loss_method() == 2) {
        y_loss = sample_rapidity_loss_from_parametrization(y_init);
    }
    return(y_loss);
}

real Glauber::sample_rapidity_loss_from_the_LEXUS_model(real y_init) const {
    const real shape_coeff = 1.0;
    real sinh_y_lrf = sinh(shape_coeff*y_init);
    real arcsinh_factor = (ran_gen_ptr.lock()->rand_uniform()
                           *(sinh(2.*shape_coeff*y_init) - sinh_y_lrf)
                           + sinh_y_lrf);
    real y_loss = 2.*y_init - asinh(arcsinh_factor)/shape_coeff;
    return(y_loss);
}

real Glauber::sample_rapidity_loss_from_parametrization(real y_init) const {
    const real slope = 0.55;
    const real upper = 0.9;
    const real transt = 2.1;
    auto y_loss = (slope*y_init
            + 0.5*(tanh(y_init - transt) + 1.)*(upper - 0.65*slope*y_init)
            - upper/2.*(1. - tanh(transt)));
    return(y_loss);
}

// sample y from exp[(y - (0.5 (yt + yp)))/2]/(4.` Sinh[0.25` yp - 0.25` yt]),
// the new rapidity of the baryon number from the right moving particle
// after the collision in the lab frame
real Glauber::sample_junction_rapidity_right(real y_left, real y_right) const {
    real y = -2.*(-0.25*y_right - 0.25*y_left
                  - 1.*log(2.*(ran_gen_ptr.lock()->rand_uniform()
                               + 0.5*exp(-0.25*y_right + 0.25*y_left)
                                 /sinh(0.25*y_right - 0.25*y_left))
                           *sinh(0.25*y_right - 0.25*y_left)));
    return(y);
}

// sample y from exp[-(y - (0.5 (yt + yp)))/2]/(4.` Sinh[0.25` yp - 0.25` yt]),
// the new rapidity of the baryon number from the left moving particle
// after the collision in the lab frame
real Glauber::sample_junction_rapidity_left(real y_left, real y_right) const {
    real y = 2.*(0.25*y_right + 0.25*y_left
                 - log(2.*(ran_gen_ptr.lock()->rand_uniform()
                           + 0.5*exp(-0.25*y_right + 0.25*y_left)
                             /sinh(0.25*y_right - 0.25*y_left))
                       *sinh( 0.25 * y_right - 0.25 * y_left)));
    return(y);
}

}
