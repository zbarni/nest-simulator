/*
 *  ShouvalNeuron.cpp
 *
 *  This file is part of NEST.
 *
 *  Copyright (C) 2004 The NEST Initiative
 *
 *  NEST is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  NEST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with NEST.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include "shouval_neuron.h"

#ifdef HAVE_GSL

// C++ includes:
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <limits>

// Includes from libnestutil:
#include "dict_util.h"
#include "numerics.h"

// Includes from nestkernel:
#include "event.h"
#include "exceptions.h"
#include "kernel_manager.h"
#include "universal_data_logger_impl.h"

// Includes from sli:
#include "dict.h"
#include "dictutils.h"
#include "doubledatum.h"
#include "integerdatum.h"

/* ----------------------------------------------------------------
 * Recordables map
 * ---------------------------------------------------------------- */

nest::RecordablesMap< nest::ShouvalNeuron > nest::ShouvalNeuron::recordablesMap_;

namespace nest // template specialization must be placed in namespace
{
// Override the create() method with one call to RecordablesMap::insert_()
// for each quantity to be recorded.
template <>
void
RecordablesMap< ShouvalNeuron >::create()
{
  // use standard names whereever you can for consistency!
  insert_( names::V_m, &ShouvalNeuron::get_y_elem_< ShouvalNeuron::State_::V_M > );  // TODO change this to V_m
  insert_( names::g_ex, &ShouvalNeuron::get_y_elem_< ShouvalNeuron::State_::G_EXC > );
  insert_( names::g_in, &ShouvalNeuron::get_y_elem_< ShouvalNeuron::State_::G_INH > );

  // TODO check if really needed?
//  insert_("g_ex__X__spikeExcRec0", &ShouvalNeuron::get_g_ex__X__spikeExcRec0);
//  insert_("g_ex__X__spikeExcRec1", &ShouvalNeuron::get_g_ex__X__spikeExcRec1);
//  insert_("g_in__X__spikeInh", &ShouvalNeuron::get_g_in__X__spikeInh);
}
}

extern "C" inline int
nest::ShouvalNeuron_dynamics( double, const double y[], double f[], void* pnode )
{
  // a shorthand
  typedef nest::ShouvalNeuron::State_ S;

  // get access to node so we can almost work as in a member function
  assert( pnode );
  const nest::ShouvalNeuron& node = *( reinterpret_cast< nest::ShouvalNeuron* >( pnode ) );

  const bool is_refractory = node.S_.r_ > 0;

  // y[] here is---and must be---the state vector supplied by the integrator,
  // not the state vector in the node, node.S_.y[].

  // The following code is verbose for the sake of clarity. We assume that a
  // good compiler will optimize the verbosity away ...

  // Clamp membrane potential to V_reset while refractory, otherwise bound
  // it to V_th.
  const double V = is_refractory ? node.P_.V_reset_ : std::min( y[ S::V_M ], node.P_.V_th_ );

  const double I_syn_exc = y[ S::G_EXC ] * ( V - node.P_.E_ex );
  const double I_syn_inh = y[ S::G_INH ] * ( V - node.P_.E_in );
  const double I_L = node.P_.g_L * ( V - node.P_.E_L );

  // V dot
  f[ 0 ] = is_refractory ? 0.0 : ( -I_L + node.B_.I_stim_ + node.P_.I_e - I_syn_exc - I_syn_inh ) / node.P_.C_m;

  f[ 1 ] = -y[ S::G_EXC ] / node.P_.tau_synE;
  f[ 2 ] = -y[ S::G_INH ] / node.P_.tau_synI;

  // TODO change here
//  // ode_state[] here is---and must be---the state vector supplied by the integrator,
//  // not the state vector in the node, node.S_.ode_state[].
//
//  const double I_syn_exc_rec1 = ode_state[ State_::g_ex__X__spikeExcRec1 ] * (-ode_state[ State_::V_m ]
//                                  + node.P_.E_ex );
//  const double I_syn_exc_rec0 = ode_state[ State_::g_ex__X__spikeExcRec0 ] * (-ode_state[ State_::V_m ]
//                                  + node.P_.E_ex );
//  const double I_syn_inh = ode_state[ State_::g_in__X__spikeInh ] * ( -ode_state[ State_::V_m ] + node.P_.E_in );
//  const double I_L = node.P_.g_L * ( -ode_state[ State_::V_m ] + node.P_.E_L );
//
//  f[State_::V_m] = (node.get_I_e() + node.B_.I_stim_grid_sum_ + I_L + I_syn_exc_rec0 + I_syn_inh + I_syn_exc_rec1)
//    / node.get_C_m();
//
//  f[State_::g_ex__X__spikeExcRec0] = -(ode_state[State_::g_ex__X__spikeExcRec0]) / node.get_tau_syn_ex();
//  f[State_::g_ex__X__spikeExcRec1] = -(ode_state[State_::g_ex__X__spikeExcRec1]) / node.get_tau_syn_ex_rec1();
//  f[State_::g_in__X__spikeInh] = -(ode_state[State_::g_in__X__spikeInh]) / node.get_tau_syn_in();

  return GSL_SUCCESS;
}

/* ----------------------------------------------------------------
 * Default constructors defining default parameters and state
 * ---------------------------------------------------------------- */

nest::ShouvalNeuron::Parameters_::Parameters_()
  : V_th_( -55.0 )    // mV
  , V_reset_( -60.0 ) // mV
  , t_ref_( 2.0 )     // ms
  , g_L( 10. )    // nS
  , C_m( 200.0 )      // pF
  , E_ex( -5.0 )       // mV
  , E_in( -70.0 )     // mV
  , E_L( -60.0 )      // mV
  , tau_synE( 0.2 )   // ms
  , tau_synI( 2.0 )   // ms
  , I_e( 0.0 )        // pA

  , tau_syn_ex_rec0 ( 80. )
  , tau_syn_ex_rec1 ( 10. )
  , tau_syn_in ( 10. )
  , tau_w ( 40. )
  , rho ( 1 / 7. )
{
}

nest::ShouvalNeuron::State_::State_( const Parameters_& p )
  : r_( 0 )
{
  y_[ V_M ] = p.E_L;
  y_[ G_EXC ] = y_[ G_INH ] = 0;

  y_[ V_m ] = p.E_L;
  y_[ g_ex__X__spikeExcRec0 ] = 0;
  y_[ g_ex__X__spikeExcRec1 ] = 0;
  y_[ g_in__X__spikeInh ] = 0;
}

nest::ShouvalNeuron::State_::State_( const State_& s )
  : r_( s.r_ )
{
  for ( size_t i = 0; i < STATE_VEC_SIZE; ++i )
  {
    y_[ i ] = s.y_[ i ];
  }
}

nest::ShouvalNeuron::State_&
nest::ShouvalNeuron::State_::operator=( const State_& s )
{
  r_ = s.r_;
  for ( size_t i = 0; i < STATE_VEC_SIZE; ++i )
  {
    y_[ i ] = s.y_[ i ];
  }
  return *this;
}

nest::ShouvalNeuron::Variables_::Variables_( )
{
  rate = 0.;
  t_end_last_trial = 1e10;
}

nest::ShouvalNeuron::Variables_::Variables_( const Variables_& __n )
{
  RefractoryCounts = __n.RefractoryCounts;
  __h = __n.__h;
  rate = __n.rate;
  t_end_last_trial = __n.t_end_last_trial;

  __P__g_ex__X__spikeRec0__g_Rec0__X__spikeRec0 = __n.__P__g_ex__X__spikeRec0__g_Rec0__X__spikeRec0;
  __P__g_ex__X__spikeRec1__g_Rec1__X__spikeRec1 = __n.__P__g_ex__X__spikeRec1__g_Rec1__X__spikeRec1;
  __P__g_in__X__spikeInh__g_in__X__spikeInh = __n.__P__g_in__X__spikeInh__g_in__X__spikeInh;
  __rate_kernel = __n.__rate_kernel;
}
/* ----------------------------------------------------------------
 * Parameter and state extractions and manipulation functions
 * ---------------------------------------------------------------- */

void
nest::ShouvalNeuron::Parameters_::get( DictionaryDatum& d ) const
{
  def< double >( d, names::V_th, V_th_ );
  def< double >( d, names::V_reset, V_reset_ );
  def< double >( d, names::t_ref, t_ref_ );
  def< double >( d, names::g_L, g_L );
  def< double >( d, names::E_L, E_L );
  def< double >( d, names::E_ex, E_ex );
  def< double >( d, names::E_in, E_in );
  def< double >( d, names::C_m, C_m );
  def< double >( d, names::tau_syn_ex, tau_synE );
  def< double >( d, names::tau_syn_in, tau_synI );
  def< double >( d, names::I_e, I_e );
}

void
nest::ShouvalNeuron::Parameters_::set( const DictionaryDatum& d, Node* node )
{
  // allow setting the membrane potential
  updateValueParam< double >( d, names::V_th, V_th_, node );
  updateValueParam< double >( d, names::V_reset, V_reset_, node );
  updateValueParam< double >( d, names::t_ref, t_ref_, node );
  updateValueParam< double >( d, names::E_L, E_L, node );

  updateValueParam< double >( d, names::E_ex, E_ex, node );
  updateValueParam< double >( d, names::E_in, E_in, node );

  updateValueParam< double >( d, names::C_m, C_m, node );
  updateValueParam< double >( d, names::g_L, g_L, node );

  updateValueParam< double >( d, names::tau_syn_ex, tau_synE, node );
  updateValueParam< double >( d, names::tau_syn_in, tau_synI, node );

  updateValueParam< double >( d, names::I_e, I_e, node );
  if ( V_reset_ >= V_th_ )
  {
    throw BadProperty( "Reset potential must be smaller than threshold." );
  }
  if ( C_m <= 0 )
  {
    throw BadProperty( "Capacitance must be strictly positive." );
  }
  if ( t_ref_ < 0 )
  {
    throw BadProperty( "Refractory time cannot be negative." );
  }
  if ( tau_synE <= 0 || tau_synI <= 0 )
  {
    throw BadProperty( "All time constants must be strictly positive." );
  }
}

void
nest::ShouvalNeuron::State_::get( DictionaryDatum& d ) const
{
  def< double >( d, names::V_m, y_[ V_M ] ); // Membrane potential
  def< double >( d, names::g_ex, y_[ G_EXC ] );
  def< double >( d, names::g_in, y_[ G_INH ] );
}

void
nest::ShouvalNeuron::State_::set( const DictionaryDatum& d, const Parameters_&, Node* node )
{
  updateValueParam< double >( d, names::V_m, y_[ V_M ], node );
  updateValueParam< double >( d, names::g_ex, y_[ G_EXC ], node );
  updateValueParam< double >( d, names::g_in, y_[ G_INH ], node );
}

nest::ShouvalNeuron::Buffers_::Buffers_( ShouvalNeuron& n )
  : logger_( n )
  , s_( 0 )
  , c_( 0 )
  , e_( 0 )
{
  // Initialization of the remaining members is deferred to
  // init_buffers_().
}

nest::ShouvalNeuron::Buffers_::Buffers_( const Buffers_&, ShouvalNeuron& n )
  : logger_( n )
  , s_( 0 )
  , c_( 0 )
  , e_( 0 )
{
  // Initialization of the remaining members is deferred to
  // init_buffers_().
}

/* ----------------------------------------------------------------
 * Default and copy constructor for node, and destructor
 * ---------------------------------------------------------------- */

nest::ShouvalNeuron::ShouvalNeuron()
  : Shouval_Archiving_Node()
  , P_()
  , S_( P_ )
  , V_( )
  , B_( *this )
{
  recordablesMap_.create();
}

nest::ShouvalNeuron::ShouvalNeuron( const ShouvalNeuron& n )
  : Shouval_Archiving_Node( n )
  , P_( n.P_ )
  , S_( n.S_ )
  , V_( n.V_ )
  , B_( n.B_, *this )
{
}

nest::ShouvalNeuron::~ShouvalNeuron()
{
  // GSL structs may not have been allocated, so we need to protect destruction
  if ( B_.s_ )
  {
    gsl_odeiv_step_free( B_.s_ );
  }
  if ( B_.c_ )
  {
    gsl_odeiv_control_free( B_.c_ );
  }
  if ( B_.e_ )
  {
    gsl_odeiv_evolve_free( B_.e_ );
  }
}

/* ----------------------------------------------------------------
 * Node initialization functions
 * ---------------------------------------------------------------- */
//  TODO check if needed
//void ShouvalNeuron::init_state_(const Node& proto){
//  const ShouvalNeuron& pr = downcast<ShouvalNeuron>(proto);
//  S_ = pr.S_;
//}

void
nest::ShouvalNeuron::init_buffers_()
{
  B_.spike_exc_.clear(); // includes resize
  B_.spike_inh_.clear(); // includes resize
  B_.currents_.clear();  // includes resize
  Shouval_Archiving_Node::clear_history();

  B_.logger_.reset();

  B_.step_ = Time::get_resolution().get_ms();
  B_.IntegrationStep_ = B_.step_;

  if ( B_.s_ == 0 )
  {
    B_.s_ = gsl_odeiv_step_alloc( gsl_odeiv_step_rkf45, State_::STATE_VEC_SIZE );
  }
  else
  {
    gsl_odeiv_step_reset( B_.s_ );
  }

  if ( B_.c_ == 0 )
  {
    B_.c_ = gsl_odeiv_control_y_new( 1e-3, 0.0 );
  }
  else
  {
    gsl_odeiv_control_init( B_.c_, 1e-3, 0.0, 1.0, 0.0 );
  }

  if ( B_.e_ == 0 )
  {
    B_.e_ = gsl_odeiv_evolve_alloc( State_::STATE_VEC_SIZE );
  }
  else
  {
    gsl_odeiv_evolve_reset( B_.e_ );
  }

  B_.sys_.function = ShouvalNeuron_dynamics;
  B_.sys_.jacobian = NULL;
  B_.sys_.dimension = State_::STATE_VEC_SIZE;
  B_.sys_.params = reinterpret_cast< void* >( this );

  B_.I_stim_ = 0.0;
}

void
nest::ShouvalNeuron::pre_run_hook()
{
  // ensures initialization in case mm connected after Simulate
  B_.logger_.init();

  V_.RefractoryCounts_ = Time( Time::ms( P_.t_ref_ ) ).get_steps();
  // since t_ref_ >= 0, this can only fail in error
  assert( V_.RefractoryCounts_ >= 0 );

  V_.__h =nest::Time::get_resolution().get_ms();
  V_.__P__g_ex__X__spikeRec0__g_Rec0__X__spikeRec0 =std::exp((-V_.__h) / P_.tau_syn_ex_rec0);
  V_.__P__g_ex__X__spikeRec1__g_Rec1__X__spikeRec1 =std::exp((-V_.__h) / P_.tau_syn_ex_rec1);
  V_.__P__g_in__X__spikeInh__g_in__X__spikeInh =std::exp((-V_.__h) / P_.tau_syn_in);
  V_.__rate_kernel =std::exp((-V_.__h) / P_.tau_w);
}

/* ----------------------------------------------------------------
 * Update and spike handling functions
 * ---------------------------------------------------------------- */

void nest::ShouvalNeuron::evolve_synaptic_activation_traces( nest::Time const & origin, const long lag )
{
  std::ostringstream msg;

  // iterate over exc and inh
  for(auto& it : V_.activeSources)
  {
    // branch on E/I
    if (V_.preSynWeights[it] >= 0.) {
      // branch on receptor type 0/1
      if (V_.isExcRec1[it]) {
        V_.preSynActivationTraces[it] *= V_.__P__g_ex__X__spikeRec1__g_Rec1__X__spikeRec1;
      }
      else {
        V_.preSynActivationTraces[it] *= V_.__P__g_ex__X__spikeRec0__g_Rec0__X__spikeRec0;
      }
    }
    else {
      V_.preSynActivationTraces[it] *= V_.__P__g_in__X__spikeInh__g_in__X__spikeInh;
    }
  }

  // evolve synaptic activations s_ij, i.e., process any spikes and add corresponding scaled delta impulses
  msg << "[update] iterating over spike events... \n";
  for(auto it = V_.spikeEvents.begin(); it != V_.spikeEvents.end();)
  {
    msg << "[update] spike from source " << it->id_ << ", scheduled @ " << it->deliveryTime
        << "; origin: " << origin.get_steps() << "; lag: " << lag << "\n";
    if (it->deliveryTime == origin.get_steps() + lag - 1)
    {
      //            msg << "[update] Adding delta impulse from source " << it->id_ << " => "
      //                << V_.preSynActivationTraces[it->id_] << "\n";

      // spike from excitatory neuron
      if (it->w >= 0.) {
        if (B_.spikeExc_grid_sum_ == 0) {
          // we allow for a 1 step difference
          if (it->one_step_mercy)
          {
            //                        LOG(nest::M_DEBUG, "ShouvalNeuron::update()", msg.str());
            LOG(nest::M_ERROR, "ShouvalNeuron:update():", "B_.spikeExc_grid_sum_ == 0");
            assert(B_.spikeExc_grid_sum_ != 0);
          }
          else
          {
            it->one_step_mercy = true;
            ++it;
            msg << "[update] compensating for 1 step difference in spike delivery...\n";
          }
        }
        else
        {
          V_.preSynActivationTraces[it->id_] += P_.rho * (1 - V_.preSynActivationTraces[it->id_]);
          it = V_.spikeEvents.erase(it);
        }
      }
      // spike from inhibitory neuron
      else {
        if (B_.spikeInh_grid_sum_ == 0)
        {
          // we allow for a 1 step difference
          if (it->one_step_mercy)
          {
            //                        LOG(nest::M_DEBUG, "ShouvalNeuron::update()", msg.str());
            LOG(nest::M_ERROR, "ShouvalNeuron:update():", "B_.spikeInh_grid_sum_ == 0");
            assert(B_.spikeInh_grid_sum_ != 0);
          }
          else
          {
            it->one_step_mercy = true;
            ++it;
            msg << "[update] compensating for 1 step difference in spike delivery...\n";
          }
        }
        else
        {
          V_.preSynActivationTraces[it->id_] += P_.rho * (1 - V_.preSynActivationTraces[it->id_]);
          it = V_.spikeEvents.erase(it);
        }
      }
    }
    else
    {
      ++it;
    }
  }

  // compute g_ij = sum_j { W_ij * s_j }
  double tmp_w_s_exc_rec0 = 0;
  double tmp_w_s_exc_rec1 = 0;
  double tmp_w_s_inh = 0;
  for(auto& it : V_.activeSources)
  {
    if (V_.preSynWeights[it] >= 0.)
    {
      if (V_.isExcRec1[it])
      {
        tmp_w_s_exc_rec1 += V_.preSynWeights[it] * V_.preSynActivationTraces[it];
      }
      else
      {
        tmp_w_s_exc_rec0 += V_.preSynWeights[it] * V_.preSynActivationTraces[it];
      }
    }
    else
    {
      tmp_w_s_inh += V_.preSynWeights[it] * V_.preSynActivationTraces[it];
    }
  }

  // update conductance after summing over the syn activations * w
  S_.y_[State_::g_ex__X__spikeExcRec0] = tmp_w_s_exc_rec0;
  S_.y_[State_::g_ex__X__spikeExcRec1] = tmp_w_s_exc_rec1;
  S_.y_[State_::g_in__X__spikeInh] = -tmp_w_s_inh;  // ensure the conductance is still positive!

}


void
nest::ShouvalNeuron::update( Time const& origin, const long from, const long to )
{

  assert( to >= 0 && ( delay ) from < kernel().connection_manager.get_min_delay() );
  assert( from < to );

  for ( long lag = from; lag < to; ++lag )
  {

    double t = 0.0;

    // numerical integration with adaptive step size control:
    // ------------------------------------------------------
    // gsl_odeiv_evolve_apply performs only a single numerical
    // integration step, starting from t and bounded by step;
    // the while-loop ensures integration over the whole simulation
    // step (0, step] if more than one integration step is needed due
    // to a small integration step size;
    // note that (t+IntegrationStep > step) leads to integration over
    // (t, step] and afterwards setting t to step, but it does not
    // enforce setting IntegrationStep to step-t; this is of advantage
    // for a consistent and efficient integration across subsequent
    // simulation intervals
    while ( t < B_.step_ )
    {
      const int status = gsl_odeiv_evolve_apply( B_.e_,
        B_.c_,
        B_.s_,
        &B_.sys_,             // system of ODE
        &t,                   // from t
        B_.step_,             // to t <= step
        &B_.IntegrationStep_, // integration step size
        S_.y_ );              // neuronal state
      if ( status != GSL_SUCCESS )
      {
        throw GSLSolverFailure( get_name(), status );
      }
    }

    S_.y_[ State_::G_EXC ] += B_.spike_exc_.get_value( lag );
    S_.y_[ State_::G_INH ] += B_.spike_inh_.get_value( lag );

    // absolute refractory period
    if ( S_.r_ )
    { // neuron is absolute refractory
      --S_.r_;
      S_.y_[ State_::V_M ] = P_.V_reset_;
    }
    else
      // neuron is not absolute refractory
      if ( S_.y_[ State_::V_M ] >= P_.V_th_ )
      {
        S_.r_ = V_.RefractoryCounts_;
        S_.y_[ State_::V_M ] = P_.V_reset_;

        set_spiketime( Time::step( origin.get_steps() + lag + 1 ) );

        SpikeEvent se;
        se.set_sender( *this );
        se.set_sender_node_id( this->get_node_id() );
        kernel().event_delivery_manager.send( *this, se, lag );
      }

    // set new input current
    B_.I_stim_ = B_.currents_.get_value( lag );

    // log state data
    B_.logger_.record_data( origin.get_steps() + lag );
  }
}

void
nest::ShouvalNeuron::handle( SpikeEvent& e )
{
  assert( e.get_delay_steps() > 0 );

  std::cout << "Spike handling from sender ID: " << e.is_valid() << " 1  " << std::endl << std::flush << std::flush;
  std::cout << "Spike handling from sender ID: " << e.receiver_is_valid() << " 2  " << std::endl << std::flush << std::flush;
  std::cout << "Spike handling from sender ID: " << e.sender_is_valid() << " 3  " << std::endl << std::flush << std::flush;

  if ( e.get_weight() > 0.0 )
  {
    B_.spike_exc_.add_value( e.get_rel_delivery_steps( kernel().simulation_manager.get_slice_origin() ),
      e.get_weight() * e.get_multiplicity() );
  }
  else
  {
    B_.spike_inh_.add_value( e.get_rel_delivery_steps( kernel().simulation_manager.get_slice_origin() ),
      -e.get_weight() * e.get_multiplicity() );
  }

  std::cout << " e.get_sender_node_id(): " <<  e.get_sender_node_id() << std::endl << std::flush;

  return;

  // NEW version

  const double weight = e.get_weight();
  const double multiplicity = e.get_multiplicity();
  const nest::rport port = e.get_rport();

  // ignore every spike from previous trial @critical
  if (e.get_stamp().get_steps() + e.get_delay_steps() >=
    kernel().simulation_manager.get_clock().delay_ms_to_steps(V_.t_end_last_trial))
  {
    return;
  }

  if (weight < 0.0) { // inhibitory
    if ( port == 0 ) {
      throw BadProperty("Because there are only I -> M connections with tau_syn_MI = 10ms = tau_syn_ex_rec1,"
        "inhibitory connections must be made onto receptor 1 (for now).");
    }
    // this includes the delay
    long deliveryTime = e.get_rel_delivery_steps(kernel().simulation_manager.get_slice_origin());

    TraceTracker_ spikeEventStruct = {e.get_stamp().get_steps() + e.get_delay_steps() - 2,
      weight, e.get_sender_node_id(), port, false, true};

    V_.spikeEvents.push_back(spikeEventStruct);
    V_.activeSources.insert(e.get_sender_node_id());
    V_.preSynWeights[e.get_sender_node_id()] = weight;
//    get_spikeInh().add_value(deliveryTime, -1 * weight * multiplicity);
    B_.spike_inh_.add_value(deliveryTime, -1 * weight * multiplicity);
  }

  if (weight >= 0.0) { // excitatory
    // this includes the delay
    long deliveryTime = e.get_rel_delivery_steps(kernel().simulation_manager.get_slice_origin());

    std::ostringstream msg;
    msg << "\n\t\tweight: " << weight
        << "\n\t\t sender id: " << e.get_sender_node_id()
        << "\n\t\t TIME: (event tstamp) " << e.get_stamp().get_steps()
//        << "\n\t\t slice_origin(): " << kernel().simulation_manager.get_slice_origin()
        << "\n\t\t get_delay_steps(): " << e.get_delay_steps()
        << "\n\t\t deliveryTime: " << deliveryTime
        << "\n\t\t absolut deliveryTime (! precise): " << e.get_stamp().get_steps() + e.get_delay_steps() - 2
        << "\n\t\t isLGN (port 1, tau_syn_rec1)? : " << (bool) (port > 0)
        << "\n";

    TraceTracker_ spikeEventStruct = {e.get_stamp().get_steps() + e.get_delay_steps() - 2,
      weight, e.get_sender_node_id(), port, false, true};

    if ( port > 0 ) {
      V_.isExcRec1[e.get_sender_node_id()] = true;
    }

    V_.spikeEvents.push_back(spikeEventStruct);
    V_.preSynWeights[e.get_sender_node_id()] = weight;
    V_.activeSources.insert(e.get_sender_node_id());
//    get_spikeExc().add_value(deliveryTime, weight * multiplicity);
    B_.spike_exc_.add_value(deliveryTime, weight * multiplicity);
  }
}

void
nest::ShouvalNeuron::handle( CurrentEvent& e )
{
  assert( e.get_delay_steps() > 0 );

  const double c = e.get_current();
  const double w = e.get_weight();

  B_.currents_.add_value( e.get_rel_delivery_steps( kernel().simulation_manager.get_slice_origin() ), w * c );
}

void
nest::ShouvalNeuron::handle( DataLoggingRequest& e )
{
  B_.logger_.handle( e );
}

#endif // HAVE_GSL
