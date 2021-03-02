/*
 *  lin_rate.h
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

#ifndef LIN_RATE_H
#define LIN_RATE_H

// Includes from models:
#include "rate_neuron_ipn.h"
#include "rate_neuron_ipn_impl.h"
#include "rate_neuron_opn.h"
#include "rate_neuron_opn_impl.h"
#include "rate_transformer_node.h"
#include "rate_transformer_node_impl.h"

namespace nest
{

/** @BeginDocumentation
@ingroup Neurons
@ingroup rate

Name: lin_rate - Linear rate model

Description:

lin_rate is an implementation of a linear rate model with
input function \f$ input(h) = g * h \f$.
The model supports multiplicative coupling which can
be switched on and off via the boolean parameter mult_coupling
(default=false). In case multiplicative coupling is actived
the excitatory input of the model is multiplied with the function
\f$ mult\_coupling\_ex(rate) = g_{ex} * ( \theta_{ex} - rate ) \f$
and the inhibitory input is multiplied with the function
\f$ mult\_coupling\_in(rate) = g_{in} * ( \theta_{in} + rate ) \f$.

The model supports connections to other rate models with either zero or
non-zero delay, and uses the secondary_event concept introduced with
the gap-junction framework.

Parameters:

The following parameters can be set in the status dictionary.
\verbatim embed:rst
===============  ======= ==================================================
 rate            real    Rate (unitless)
 tau             ms      Time constant of rate dynamics
 lambda          real    Passive decay rate
 mu              real    Mean input
 sigma           real    Noise parameter
 g               real    Gain parameter
 mult_coupling   boolean Switch to enable/disable multiplicative coupling
 g_ex            real    Linear factor in multiplicative coupling
 g_in            real    Linear factor in multiplicative coupling
 theta_ex        real    Shift in multiplicative coupling
 theta_in        real    Shift in multiplicative coupling
 rectify_output  boolean Switch to restrict rate to values >= 0
===============  ======= ==================================================
\endverbatim

References:

\verbatim embed:rst
.. [1] Hahne J, Dahmen D, Schuecker J, Frommer A, Bolten M, Helias M, Diesmann
       M (2017). Integration of continuous-time dynamics in a spiking neural
       network simulator. Frontiers in Neuroinformatics, 11:34.
       DOI: https://doi.org/10.3389/fninf.2017.00034
.. [2] Hahne J, Helias M, Kunkel S, Igarashi J, Bolten M, Frommer A, Diesmann M
       (2015). A unified framework for spiking and gap-junction interactions
       in distributed neuronal network simulations.
       Frontiers Neuroinformatics, 9:22.
       DOI: https://doi.org/10.3389/fninf.2015.00022
\endverbatim

Sends: InstantaneousRateConnectionEvent, DelayedRateConnectionEvent

Receives: InstantaneousRateConnectionEvent, DelayedRateConnectionEvent,
DataLoggingRequest

Author: David Dahmen, Jan Hahne, Jannis Schuecker

SeeAlso: rate_connection_instantaneous, rate_connection_delayed
*/

class nonlinearities_lin_rate
{
private:
  /** gain factor of gain function */
  double g_;
  /** linear factor in multiplicative excitatory coupling*/
  double g_ex_;
  /** linear factor in multiplicative inhibitory coupling*/
  double g_in_;
  /** offset in multiplicative coupling*/
  double theta_ex_;
  double theta_in_;

public:
  /** sets default parameters */
  nonlinearities_lin_rate()
    : g_( 1.0 )
    , g_ex_( 1.0 )
    , g_in_( 1.0 )
    , theta_ex_( 0.0 )
    , theta_in_( 0.0 )
  {
  }

  void get( DictionaryDatum& ) const; //!< Store current values in dictionary
  void set( const DictionaryDatum& ); //!< Set values from dicitonary

  double input( double h );               // non-linearity on input
  double mult_coupling_ex( double rate ); // factor of multiplicative coupling
  double mult_coupling_in( double rate ); // factor of multiplicative coupling
};

inline double
nonlinearities_lin_rate::input( double h )
{
  return g_ * h;
}

inline double
nonlinearities_lin_rate::mult_coupling_ex( double rate )
{
  return g_ex_ * ( theta_ex_ - rate );
}

inline double
nonlinearities_lin_rate::mult_coupling_in( double rate )
{
  return g_in_ * ( theta_in_ + rate );
}

typedef rate_neuron_ipn< nest::nonlinearities_lin_rate > lin_rate_ipn;
typedef rate_neuron_opn< nest::nonlinearities_lin_rate > lin_rate_opn;
typedef rate_transformer_node< nest::nonlinearities_lin_rate > rate_transformer_lin;

template <>
void RecordablesMap< lin_rate_ipn >::create();
template <>
void RecordablesMap< lin_rate_opn >::create();
template <>
void RecordablesMap< rate_transformer_lin >::create();


} // namespace nest


#endif /* #ifndef LIN_RATE_H */
