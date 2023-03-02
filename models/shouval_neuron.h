#ifndef SHOUVAL_NEURON_H
#define SHOUVAL_NEURON_H

// Generated includes:
#include "config.h"

#ifdef HAVE_GSL

// External includes:
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

// Includes from nestkernel:
//#include "archiving_node.h"
#include "shouval_archiving_node.h"
#include "connection.h"
#include "event.h"
#include "nest_types.h"
#include "recordables_map.h"
#include "ring_buffer.h"
#include "universal_data_logger.h"

namespace nest
{
/**
 * Function computing right-hand side of ODE for GSL solver.
 * @note Must be declared here so we can befriend it in class.
 * @note Must have C-linkage for passing to GSL. Internally, it is
 *       a first-class C++ function, but cannot be a member function
 *       because of the C-linkage.
 * @note No point in declaring it inline, since it is called
 *       through a function pointer.
 * @param void* Pointer to model neuron instance.
 */
extern "C" int ShouvalNeuron_dynamics( double, const double*, double*, void* );


//class ShouvalNeuron : public ArchivingNode
class ShouvalNeuron : public Shouval_Archiving_Node
{

public:
  ShouvalNeuron();
  ShouvalNeuron( const ShouvalNeuron& );
  ~ShouvalNeuron();

  /**
   * Import sets of overloaded virtual functions.
   * @see Technical Issues / Virtual Functions: Overriding, Overloading, and
   * Hiding
   */
  using Node::handle;
  using Node::handles_test_event;

  port send_test_event( Node&, rport, synindex, bool ) override;

  void handle( SpikeEvent& ) override;
  void handle( CurrentEvent& ) override;
  void handle( DataLoggingRequest& ) override;

  port handles_test_event( SpikeEvent&, rport ) override;
  port handles_test_event( CurrentEvent&, rport ) override;
  port handles_test_event( DataLoggingRequest&, rport ) override;

  void get_status( DictionaryDatum& ) const override;
  void set_status( const DictionaryDatum& ) override;

private:
  void init_buffers_() override;
  void pre_run_hook() override;
  void update( Time const&, const long, const long ) override;
  void evolve_synaptic_activation_traces( nest::Time const &, const long );

  // END Boilerplate function declarations ----------------------------

  // Friends --------------------------------------------------------

  // make dynamics function quasi-member
  friend int ShouvalNeuron_dynamics( double, const double*, double*, void* );

  // The next two classes need to be friends to access the State_ class/member
  friend class RecordablesMap< ShouvalNeuron >;
  friend class UniversalDataLogger< ShouvalNeuron >;


  struct TraceTracker_{
    long deliveryTime;
    double w;
    unsigned long id_;
    nest::rport port;
    // this allows a one simulation step difference in spike delivery, in order to allow correct processing;
    // spike_generator and poisson_generators call the handle() function with a 1 step difference, and we must
    // compensate for this by waiting for 1 extra step for the spike arrival during the update().
    bool one_step_mercy;
    bool is_exc;

  };

private:
  // ----------------------------------------------------------------

  //! Model parameters
  struct Parameters_
  {
    double V_th_;    //!< Threshold Potential in mV
    double V_reset_; //!< Reset Potential in mV
    double t_ref_;   //!< Refractory period in ms
//    double g_L;      //!< Leak Conductance in nS
//    double C_m;      //!< Membrane Capacitance in pF
//    double E_ex;     //!< Excitatory reversal Potential in mV
//    double E_in;     //!< Inhibitory reversal Potential in mV
//    double E_L;      //!< Leak reversal Potential (aka resting potential) in mV
//    double tau_synE; //!< Time constant for excitatory synaptic kernel in ms
//    double tau_synI; //!< Time constant for inhibitory synaptic kernel in ms
//    double I_e;      //!< Constant Current in pA

    bool logOutput;

    //!  Threshold Potential
    double V_th;

    //!  Reset Potential
    double V_reset;

    //!  Refractory period
    double t_ref;

    //!  Leak Conductance, 0.001 microS
    double g_L;

    //!  Membrane Capacitance, 20*g_L in nF
    double C_m;

    //!  Excitatory reversal Potential
    double E_ex;

    //!  Inhibitory reversal Potential
    double E_in;

    //!  Leak reversal Potential (aka resting potential)
    double E_L;

    double tau_syn_ex_rec1;

    //!  Synaptic Time Constant Excitatory Synapse
    double tau_syn_ex_rec0;

    //!  Synaptic Time Constant for Inhibitory Synapse
    double tau_syn_in;

    //! Rate estimation (filtering) time constant
    double tau_w;

    //!  Fractional change of synaptic activation
    double rho;

    //!  constant external input current
    double I_e;

    double __gsl_error_tol;

    Parameters_(); //!< Sets default parameter values

    void get( DictionaryDatum& ) const;             //!< Store current values in dictionary
    void set( const DictionaryDatum&, Node* node ); //!< Set values from dicitonary
  };

public:
  // ----------------------------------------------------------------

  /**
   * State variables of the model.
   * @note Copy constructor required because of C-style array.
   */
  struct State_
  {
    //! Symbolic indices to the elements of the state vector y
    enum StateVecElems
    {
//      V_M = 0,  // TODO del
//      G_EXC,
//      G_INH,

      V_m,
      g_ex__X__spikeExcRec1,
      g_ex__X__spikeExcRec0,
      g_in__X__spikeInh,
      STATE_VEC_SIZE
    };

    //! neuron state, must be C-array for GSL solver
    double y_[ STATE_VEC_SIZE ];
//    double ode_state[STATE_VEC_SIZE];
    int r_; //!< number of refractory steps remaining

    State_( const Parameters_& ); //!< Default initialization
    State_( const State_& );

    State_& operator=( const State_& );

    void get( DictionaryDatum& ) const;
    void set( const DictionaryDatum&, const Parameters_&, Node* );
  };

  // ----------------------------------------------------------------

private:
  /**
   * Buffers of the model.
   */
  struct Buffers_
  {
    Buffers_( ShouvalNeuron& );                  //!< Sets buffer pointers to 0
    Buffers_( const Buffers_&, ShouvalNeuron& ); //!< Sets buffer pointers to 0

    //! Logger for all analog data
    UniversalDataLogger< ShouvalNeuron > logger_;

    /** buffers and sums up incoming spikes/currents */
    RingBuffer spike_exc_;
    RingBuffer spike_inh_;
    RingBuffer currents_;

    double spikeInh_grid_sum_;
    double spikeExc_grid_sum_;
    double I_stim_grid_sum_;

    /** GSL ODE stuff */
    gsl_odeiv_step* s_;    //!< stepping function
    gsl_odeiv_control* c_; //!< adaptive stepsize control function
    gsl_odeiv_evolve* e_;  //!< evolution function
    gsl_odeiv_system sys_; //!< struct describing system

    // Since IntergrationStep_ is initialized with step_, and the resolution
    // cannot change after nodes have been created, it is safe to place both
    // here.
    double step_;            //!< step size in ms
    double IntegrationStep_; //!< current integration time step, updated by GSL

    /**
     * Input current injected by CurrentEvent.
     * This variable is used to transport the current applied into the
     * _dynamics function computing the derivative of the state vector.
     * It must be a part of Buffers_, since it is initialized once before
     * the first simulation, but not modified before later Simulate calls.
     */
    double I_stim_;
  };

  // ----------------------------------------------------------------

  /**
   * Internal variables of the model.
   */
  struct Variables_
  {
    int RefractoryCounts_;  // TODO del

    // recording variables, for debugging and visualization
    std::vector<double> recorded_times;
    std::vector<long> recorded_input_ids;
    std::vector<double> recorded_activation_traces;
    std::vector<double> recorded_rates;

    double preSynActivationTraces[10000] = {};  // synaptic activation traces, for each incoming synapse
    double preSynWeights[10000] = {}; // weights of the current/last spike from each synapse
    bool isExcRec1[10000] = {};

    std::list<TraceTracker_> spikeEvents;  // TODO add comm
    std::set<unsigned long> activeSources;  // TODO add comm
                                             //        std::set<unsigned long> activeSourcesRecType1;

    double t_end_last_trial;  // end timepoint of last trial. Used to ignore all incoming spikes from prev. trial.

    //!  refractory time in steps
    long RefractoryCounts;

    //! estimated rate (filtered)
    double rate;
    double __rate_kernel;
    double __h;
    double __P__g_ex__X__spikeRec1__g_Rec1__X__spikeRec1;
    double __P__g_ex__X__spikeRec0__g_Rec0__X__spikeRec0;
    double __P__g_in__X__spikeInh__g_in__X__spikeInh;

    Variables_( ); //!< Default initialization
    Variables_( const Variables_& );

    void get( DictionaryDatum& ) const;
    void set( const DictionaryDatum&, const Parameters_&, Node* );
  };

  // Access functions for UniversalDataLogger -------------------------------

  //! Read out state vector elements, used by UniversalDataLogger
  template < State_::StateVecElems elem >
  double
  get_y_elem_() const
  {
    return S_.y_[ elem ];
  }

  // ----------------------------------------------------------------
  // TODO @check if needed
//  /* getters/setters for functions */
//  inline double get_I_syn_exc() const {
//    return S_.ode_state[State_::g_ex__X__spikeExcRec0] * (S_.ode_state[State_::V_m] - P_.E_ex);
//  }
//
//  inline double get_I_syn_inh() const {
//    return S_.ode_state[State_::g_in__X__spikeInh] * (S_.ode_state[State_::V_m] - P_.E_in);
//  }
//
//  inline double get_g_ex() const {
//    /**
//         * This is the conductance g_ex...
//     */
//    return (S_.ode_state[State_::g_ex__X__spikeExcRec0] * (S_.ode_state[State_::V_m] - P_.E_ex))
//      / (S_.ode_state[State_::V_m] - P_.E_ex);
//  }
//
//  inline double get_g_in() const {
//    /**
//         * This is the conductance g_ex...
//     */
//    return (S_.ode_state[State_::g_in__X__spikeInh] * (S_.ode_state[State_::V_m] - P_.E_in)) /
//      (S_.ode_state[State_::V_m] - P_.E_in);
//  }
//
//  inline double get_I_leak() const {
//    return P_.g_L * (S_.ode_state[State_::V_m] - P_.E_L);
//  }

  Parameters_ P_;
  State_ S_;
  Variables_ V_;
  Buffers_ B_;

  //! Mapping of recordables names to access functions
  static RecordablesMap< ShouvalNeuron > recordablesMap_;
};


inline port
nest::ShouvalNeuron::send_test_event( Node& target, rport receptor_type, synindex, bool )
{
  SpikeEvent e;
  e.set_sender( *this );
  return target.handles_test_event( e, receptor_type );
}

inline port
ShouvalNeuron::handles_test_event( SpikeEvent&, rport receptor_type )
{
//  if ( receptor_type != 0 )
//  if ( receptor_type < 0 || receptor_type > 1 )
  if ( receptor_type < 0)
  {
    throw UnknownReceptorType( receptor_type, get_name() );
  }
  return receptor_type;
}

inline port
ShouvalNeuron::handles_test_event( CurrentEvent&, rport receptor_type )
{
  if ( receptor_type != 0 )
  {
    throw UnknownReceptorType( receptor_type, get_name() );
  }
  return 0;
}

inline port
ShouvalNeuron::handles_test_event( DataLoggingRequest& dlr, rport receptor_type )
{
  if ( receptor_type != 0 )
  {
    throw UnknownReceptorType( receptor_type, get_name() );
  }
  return B_.logger_.connect_logging_device( dlr, recordablesMap_ );
}


inline void
ShouvalNeuron::get_status( DictionaryDatum& d ) const
{
  P_.get( d );
  S_.get( d );
  V_.get( d );
//  ArchivingNode::get_status( d );
  Shouval_Archiving_Node::get_status( d );

  ( *d )[ names::recordables ] = recordablesMap_.get_list();
}

inline void
ShouvalNeuron::set_status( const DictionaryDatum& d )
{
  Parameters_ ptmp = P_;     // temporary copy in case of errors
  ptmp.set( d, this );       // throws if BadProperty

  State_ stmp = S_;          // temporary copy in case of errors
  stmp.set( d, ptmp, this ); // throws if BadProperty

  Variables_ vtmp = V_;          // temporary copy in case of errors
  vtmp.set( d, ptmp, this ); // throws if BadProperty

  // We now know that (ptmp, stmp) are consistent. We do not
  // write them back to (P_, S_) before we are also sure that
  // the properties to be set in the parent class are internally
  // consistent.
//  ArchivingNode::set_status( d );
  Shouval_Archiving_Node::set_status( d );

  // if we get here, temporaries contain consistent set of properties
  P_ = ptmp;
  S_ = stmp;
  V_ = vtmp;
}

} // namespace

#endif // HAVE_GSL
#endif // SHOUVAL_NEURON_H
