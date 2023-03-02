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
  insert_( names::V_m, &ShouvalNeuron::get_y_elem_< ShouvalNeuron::State_::V_m > );  // TODO change this to V_m
//  insert_( names::g_ex, &ShouvalNeuron::get_y_elem_< ShouvalNeuron::State_::G_EXC > );
//  insert_( names::g_in, &ShouvalNeuron::get_y_elem_< ShouvalNeuron::State_::G_INH > );

  // TODO check if really needed?
  insert_("g_ex__X__spikeExcRec0", &ShouvalNeuron::get_y_elem_< ShouvalNeuron::State_::g_ex__X__spikeExcRec0 > );
  insert_("g_ex__X__spikeExcRec1", &ShouvalNeuron::get_y_elem_< ShouvalNeuron::State_::g_ex__X__spikeExcRec1 > );
  insert_("g_in__X__spikeInh", &ShouvalNeuron::get_y_elem_< ShouvalNeuron::State_::g_in__X__spikeInh > );
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

  // Clamp membrane potential to V_reset while refractory, otherwise bound
  // it to V_th.
//  const double V = is_refractory ? node.P_.V_reset_ : std::min( y[ S::V_m ], node.P_.V_th_ );
//
//  const double I_syn_exc = y[ S::G_EXC ] * ( V - node.P_.E_ex );
//  const double I_syn_inh = y[ S::G_INH ] * ( V - node.P_.E_in );
//  const double I_L = node.P_.g_L * ( V - node.P_.E_L );
//
//  // V dot
//  f[ 0 ] = is_refractory ? 0.0 : ( -I_L + node.B_.I_stim_ + node.P_.I_e - I_syn_exc - I_syn_inh ) / node.P_.C_m;
//
//  f[ 1 ] = -y[ S::G_EXC ] / node.P_.tau_synE;
//  f[ 2 ] = -y[ S::G_INH ] / node.P_.tau_synI;

  // TODO it seems that we have not been clamping the voltage to reset properly...
//  const double V = is_refractory ? node.P_.V_reset_ : std::min( y[ S::V_m ], node.P_.V_th_ );
  const double V = y[ S::V_m ];

  // this was before...
  const double I_syn_exc_rec0 = y[ S::g_ex__X__spikeExcRec0 ] * (-V + node.P_.E_ex );
  const double I_syn_exc_rec1 = y[ S::g_ex__X__spikeExcRec1 ] * (-V + node.P_.E_ex );
  const double I_syn_inh = y[ S::g_in__X__spikeInh ] * ( -V + node.P_.E_in );
  const double I_L = node.P_.g_L * ( -V + node.P_.E_L );

    // TODO it seems that we have not been clamping the voltage to reset properly...
//  f[S::V_m] = is_refractory ? 0.0 : (node.P_.I_e + node.B_.I_stim_grid_sum_ + I_L
//                                     + I_syn_exc_rec0 + I_syn_inh + I_syn_exc_rec1) / node.P_.C_m;
  f[S::V_m] = (node.P_.I_e + node.B_.I_stim_grid_sum_ + I_L + I_syn_exc_rec0 + I_syn_inh + I_syn_exc_rec1)
              / node.P_.C_m;

  f[S::g_ex__X__spikeExcRec0] = -(y[S::g_ex__X__spikeExcRec0]) / node.P_.tau_syn_ex_rec0;
  f[S::g_ex__X__spikeExcRec1] = -(y[S::g_ex__X__spikeExcRec1]) / node.P_.tau_syn_ex_rec1;
  f[S::g_in__X__spikeInh] = -(y[S::g_in__X__spikeInh]) / node.P_.tau_syn_in;

  return GSL_SUCCESS;
}

/* ----------------------------------------------------------------
 * Default constructors defining default parameters and state
 * ---------------------------------------------------------------- */

nest::ShouvalNeuron::Parameters_::Parameters_()
  : V_th( -55.0 )    // mV
  , V_reset( -61.0 ) // mV
  , t_ref( 2.0 )     // ms
  , g_L( 10. )    // nS
  , C_m( 200.0 )      // pF
  , E_ex( -5.0 )       // mV
  , E_in( -70.0 )     // mV
  , E_L( -60.0 )      // mV
  , I_e( 0.0 )        // pA

  , tau_syn_ex_rec0 ( 80. )
  , tau_syn_ex_rec1 ( 10. )
  , tau_syn_in ( 10. )
  , tau_w ( 40. )
  , rho ( 1 / 7. )
  , logOutput ( false )
{
}

nest::ShouvalNeuron::State_::State_( const Parameters_& p )
  : r_( 0 )
{
  y_[ V_m ] = p.E_L;
//  y_[ G_EXC ] = y_[ G_INH ] = 0;

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
  def< double >( d, names::V_th, V_th );
  def< double >( d, names::V_reset, V_reset );
  def< double >( d, names::t_ref, t_ref );
  def< double >( d, names::g_L, g_L );
  def< double >( d, names::E_L, E_L );
  def< double >( d, names::E_ex, E_ex );
  def< double >( d, names::E_in, E_in );
  def< double >( d, names::C_m, C_m );
  def< double >( d, names::I_e, I_e );
  
  // CS
  def<double>(d, "tau_syn_ex_rec0", tau_syn_ex_rec0);
  def<double>(d, "tau_syn_ex_rec1", tau_syn_ex_rec1);
  def<double>(d, "tau_syn_in", tau_syn_in);
  def<double>(d, "tau_w", tau_w);
  def<double>(d, "rho", rho);
  def<double>(d, "logOutput", logOutput);
}

void
nest::ShouvalNeuron::Parameters_::set( const DictionaryDatum& d, Node* node )
{
  // allow setting the membrane potential
  updateValueParam< double >( d, names::V_th, V_th, node );
  updateValueParam< double >( d, names::V_reset, V_reset, node );
  updateValueParam< double >( d, names::t_ref, t_ref, node );
  updateValueParam< double >( d, names::E_L, E_L, node );

  updateValueParam< double >( d, names::E_ex, E_ex, node );
  updateValueParam< double >( d, names::E_in, E_in, node );

  updateValueParam< double >( d, names::C_m, C_m, node );
  updateValueParam< double >( d, names::g_L, g_L, node );

  updateValueParam< double >( d, names::I_e, I_e, node );
  
  updateValueParam<double>(d, "tau_syn_ex_rec0", tau_syn_ex_rec0, node);
  updateValueParam<double>(d, "tau_syn_ex_rec1", tau_syn_ex_rec1, node);
  updateValueParam<double>(d, "tau_syn_in", tau_syn_in, node);
  updateValueParam<double>(d, "tau_w", tau_w, node);
  updateValueParam<double>(d, "rho", rho, node);
  updateValueParam<bool>(d, "logOutput", logOutput, node);
  
  if ( V_reset >= V_th )
  {
    throw BadProperty( "Reset potential must be smaller than threshold." );
  }
  if ( C_m <= 0 )
  {
    throw BadProperty( "Capacitance must be strictly positive." );
  }
  if ( t_ref < 0 )
  {
    throw BadProperty( "Refractory time cannot be negative." );
  }
}

void
nest::ShouvalNeuron::State_::get( DictionaryDatum& d ) const
{
  def< double >( d, names::V_m, y_[ V_m ] ); // Membrane potential
  def<double>(d, "g_ex__X__spikeExcRec0", y_[ g_ex__X__spikeExcRec0 ]);
  def<double>(d, "g_ex__X__spikeExcRec1", y_[ g_ex__X__spikeExcRec1 ]);
  def<double>(d, "g_in__X__spikeInh", y_[ g_in__X__spikeInh ]);
}

void
nest::ShouvalNeuron::State_::set( const DictionaryDatum& d, const Parameters_&, Node* node )
{
  updateValueParam< double >( d, names::V_m, y_[ V_m ], node );
  updateValueParam<double>(d, "g_ex__X__spikeExcRec0", y_[ g_ex__X__spikeExcRec0 ], node );
  updateValueParam<double>(d, "g_ex__X__spikeExcRec1", y_[ g_ex__X__spikeExcRec1 ], node );
  updateValueParam<double>(d, "g_in__X__spikeInh", y_[ g_in__X__spikeInh ], node );
}

void
nest::ShouvalNeuron::Variables_::get( DictionaryDatum& d ) const
{
  def< double >( d, "t_end_last_trial", t_end_last_trial );
  def< double >( d, "rate", rate );
  def<std::vector<double>>(d, "recorded_times", recorded_times);
  def<std::vector<long>>(d, "recorded_input_ids", recorded_input_ids);
  def<std::vector<double>>(d, "recorded_activation_traces", recorded_activation_traces);
  def<std::vector<double>>(d, "recorded_rates", recorded_rates);
}

void
nest::ShouvalNeuron::Variables_::set( const DictionaryDatum& d, const Parameters_& p, Node* node )
{
  updateValueParam< double >( d, "t_end_last_trial", t_end_last_trial , node );
  updateValueParam< double >( d, "rate", rate , node );
  std::vector<double> tmp_recorded_times = recorded_times;
  std::vector<long> tmp_recorded_input_ids = recorded_input_ids;
  std::vector<double> tmp_recorded_activation_traces = recorded_activation_traces;
  std::vector<double> tmp_recorded_rates = recorded_rates;

  if (p.logOutput)
  {
    tmp_recorded_times.clear();
    tmp_recorded_input_ids.clear();
    tmp_recorded_activation_traces.clear();
    tmp_recorded_rates.clear();
  }

  recorded_times = tmp_recorded_times;
  recorded_input_ids = tmp_recorded_input_ids;
  recorded_activation_traces =  tmp_recorded_activation_traces;
  recorded_rates = tmp_recorded_rates;
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

  V_.RefractoryCounts_ = Time( Time::ms( P_.t_ref ) ).get_steps();
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
//        msg << "[update] processing active source, preSynActiv Trace before any update: " << it << " => "
//            << V_.preSynActivationTraces[it] << "\n";
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
//      msg << "[update] Adding delta impulse from source " << it->id_ << " => "
//          << V_.preSynActivationTraces[it->id_] << "\n";

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
//            msg << "[update] compensating for 1 step difference in spike delivery...\n";
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
//            msg << "[update] compensating for 1 step difference in spike delivery...\n";
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
//        msg << "[update] receptor 1, presyn weight is : " << V_.preSynWeights[it]
//            << "and preSynActivationTrace is : " << V_.preSynActivationTraces[it]
//            << "\n";
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
//  std::cout << msg.str() << std::endl << std::flush;
//  std::cout << "[" << this->get_node_id() << "] @@@ updated g_ex__X__spikeExcRec = " << S_.y_[State_::g_ex__X__spikeExcRec1] << "\n"
//            << std::flush;
}


void
nest::ShouvalNeuron::update( Time const& origin, const long from, const long to )
{

  assert( to >= 0 && ( delay ) from < kernel().connection_manager.get_min_delay() );
  assert( from < to );

  std::ostringstream msg;

  for ( long lag = from; lag < to; ++lag )
  {
    double t = 0.0;

//    msg << "from " << from << " to " << to << ", lag = " << lag << "\n";

    B_.spikeInh_grid_sum_ = B_.spike_inh_.get_value(lag);
    B_.spikeExc_grid_sum_ = B_.spike_exc_.get_value(lag);
    B_.I_stim_grid_sum_ = B_.currents_.get_value(lag);


    double g_ex__X__spikeExcRec0__tmp = V_.__P__g_ex__X__spikeRec0__g_Rec0__X__spikeRec0
      * S_.y_[State_::g_ex__X__spikeExcRec0];
    double g_ex__X__spikeExcRec1__tmp = V_.__P__g_ex__X__spikeRec1__g_Rec1__X__spikeRec1
      * S_.y_[State_::g_ex__X__spikeExcRec1];
    double g_in__X__spikeInh__tmp = V_.__P__g_in__X__spikeInh__g_in__X__spikeInh
      * S_.y_[State_::g_in__X__spikeInh];
    V_.rate *= V_.__rate_kernel;

//    std::cout << "[" << this->get_node_id() << "] >>> g_ex__X__spikeExcRec = " << S_.y_[State_::g_ex__X__spikeExcRec1] << "\n"
//              << "[" << this->get_node_id() << "] >>> V_.__P__g_ex__X__spikeRec1__g_Rec1__X__spikeRec1 = " << V_.__P__g_ex__X__spikeRec1__g_Rec1__X__spikeRec1 << "\n"
//              << std::flush;
    // TODO shouldn't this be outside the loop?
    if (P_.logOutput)
    {
      V_.recorded_times.push_back( origin.get_steps() + lag );
      V_.recorded_rates.push_back( V_.rate );
    }

//    std::cout << "[" << this->get_node_id() << "] Vm before integration is " << S_.y_[State_::V_m] << std::endl << std::flush;
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

    /* replace analytically solvable variables with precisely integrated values  */
    S_.y_[State_::g_ex__X__spikeExcRec0] = g_ex__X__spikeExcRec0__tmp;
    S_.y_[State_::g_ex__X__spikeExcRec1] = g_ex__X__spikeExcRec1__tmp;
    S_.y_[State_::g_in__X__spikeInh] = g_in__X__spikeInh__tmp;

//    std::cout << "[" << this->get_node_id() << "] >>> AFTER g_ex__X__spikeExcRec = " << S_.y_[State_::g_ex__X__spikeExcRec1] << "\n" << std::flush;

    // evolves the activation traces s_i and updates the conductances scaled by the weight
    evolve_synaptic_activation_traces(origin, lag);
//    std::cout << "[" << this->get_node_id() << "] Vm after integration is " << S_.y_[State_::V_m] << std::endl << std::flush;
    // absolute refractory period
    if ( S_.r_ )
    { // neuron is absolute refractory
      --S_.r_;
      S_.y_[ State_::V_m ] = P_.V_reset;
    }
    else
      // neuron is not absolute refractory
      if ( S_.y_[ State_::V_m ] >= P_.V_th )
      {
        S_.r_ = V_.RefractoryCounts_;
        S_.y_[ State_::V_m ] = P_.V_reset;

        set_spiketime( Time::step( origin.get_steps() + lag + 1 ) );
        SpikeEvent se;
//        se.set_rport( this->get_node_id() );  // set sender port to node id so that the receiving neuron has access to it
        kernel().event_delivery_manager.send( *this, se, lag );

        //! also update the rate
        V_.rate += 1. / P_.tau_w;
      }

    // set new input current
    B_.I_stim_ = B_.currents_.get_value( lag );

    // log state data
    B_.logger_.record_data( origin.get_steps() + lag );
//    std::cout << msg.str() << std::endl << std::flush;
  }
}

void
nest::ShouvalNeuron::handle( SpikeEvent& e )
{
  assert( e.get_delay_steps() > 0 );

  const double weight = e.get_weight();
  const double multiplicity = e.get_multiplicity();
  const nest::rport encoded_rport = e.get_rport();
  nest::rport rport;
  unsigned long sender_node_id;

  // need to differentiate between spikes from other neurons and from generators
  if (encoded_rport > 1e5)
  { // from neuron
    rport = encoded_rport / (int)1e6 - 1;  // first digit encodes the rport + 1
    sender_node_id = encoded_rport % int(1e5);
  }
  else
  { // from generator
    rport = encoded_rport;
//      std::cout << "Neuron " << this->get_node_id() << " received a spike on RPORT: " << rport <<
//        "from SOURCE NEURON ? with weight: " << weight << std::endl << std::flush;
    sender_node_id = e.get_sender_node_id();  // we assume this works here, and it should hold for non-neuron nodes (?!)
  }

  assert(0 <= rport and rport <= 1);

  // ignore every spike from previous trial @critical
  if (e.get_stamp().get_steps() + e.get_delay_steps() >=
    kernel().simulation_manager.get_clock().delay_ms_to_steps(V_.t_end_last_trial))
  {
    return;
  }

//  std::cout << "Neuron " << this->get_node_id() << " received a spike on RPORT: " << rport <<
//    "from SOURCE NEURON " << sender_node_id << " with encoded PORT " << encoded_rport << std::endl << std::flush;

  // inhibitory
  if (weight < 0.0)
  {
    if ( rport == 0 )
    {
      throw BadProperty("Because there are only I -> M connections with tau_syn_MI = 10ms = tau_syn_ex_rec1,"
        "inhibitory connections must be made onto receptor 1 (for now).");
    }
    // this includes the delay
    long deliveryTime = e.get_rel_delivery_steps(kernel().simulation_manager.get_slice_origin());

    TraceTracker_ spikeEventStruct = {e.get_stamp().get_steps() + e.get_delay_steps() - 2,
      weight, sender_node_id, rport, false, true};

    V_.spikeEvents.push_back(spikeEventStruct);
    V_.activeSources.insert(sender_node_id);
    V_.preSynWeights[sender_node_id] = weight;
//    get_spikeInh().add_value(deliveryTime, -1 * weight * multiplicity);
    B_.spike_inh_.add_value(deliveryTime, -1 * weight * multiplicity);
  }

  // excitatory
  if (weight >= 0.0)
  {
    // this includes the delay
    long deliveryTime = e.get_rel_delivery_steps(kernel().simulation_manager.get_slice_origin());

    std::ostringstream msg;
    msg << "\n\t\tweight: " << weight
        << "\n\t\t sender id: " << sender_node_id
        << "\n\t\t TIME: (event tstamp) " << e.get_stamp().get_steps()
        << "\n\t\t get_delay_steps(): " << e.get_delay_steps()
        << "\n\t\t deliveryTime: " << deliveryTime
        << "\n\t\t absolut deliveryTime (! precise): " << e.get_stamp().get_steps() + e.get_delay_steps() - 2
        << "\n\t\t isLGN (port 1, tau_syn_rec1)? : " << (bool) (rport > 0)
        << "\n";

    TraceTracker_ spikeEventStruct = {e.get_stamp().get_steps() + e.get_delay_steps() - 2,
      weight, sender_node_id, rport, false, true};

    if ( rport > 0 ) {
      V_.isExcRec1[sender_node_id] = true;
    }

    V_.spikeEvents.push_back(spikeEventStruct);
    V_.preSynWeights[sender_node_id] = weight;
    V_.activeSources.insert(sender_node_id);
//    get_spikeExc().add_value(deliveryTime, weight * multiplicity);
    B_.spike_exc_.add_value(deliveryTime, weight * multiplicity);

//    std::cout << msg.str() <<  std::endl << std::flush;
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
