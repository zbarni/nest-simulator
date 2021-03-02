/*
 *  shouval_archiving_node.h
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

#ifndef SHOUVAL_ARCHIVING_NODE_H
#define SHOUVAL_ARCHIVING_NODE_H

// C++ includes:
#include <deque>

// Includes from nestkernel:
#include "histentry.h"
#include "nest_time.h"
#include "nest_types.h"
#include "archiving_node.h"
#include "synaptic_element.h"

// Includes from sli:
#include "dictdatum.h"

namespace nest
{

/**
 * \class Shouval_Archiving_Node
 * a archiving node which additionally archives parameters
 */
class Shouval_Archiving_Node : public Archiving_Node
{

public:
  /**
   * \fn Shouval_Archiving_Node()
   * Constructor.
   */
  Shouval_Archiving_Node();
  ~Shouval_Archiving_Node();

  /**
   * \fn Shouval_Archiving_Node()
   * Copy Constructor.
   */
  Shouval_Archiving_Node( const Shouval_Archiving_Node& );

  /**
   * \fn double get_LTD_value(long t)
   * Returns value in LTD history at time t
   */
//  double get_LTD_value( double t, index sender_gid ); // actually the difference now
//  double get_LTP_value( double t, index sender_gid );
//  double get_shouval_weight( double t, index sender_gid );
//  double shouval_update_weight( double t_start, double t_stop, double weight,double learn_rate, index sender_gid );
  double shouval_update_weight( double t_start, double t_stop,
                                double tau_ltp, double tau_ltd,
                                double Tp_max, double Td_max,
                                double eta_ltp, double eta_ltd,
                                double weight, double learn_rate,
                                index sender_gid );


protected:
  /**
   */
  struct connectionDataStruct
  {
    bool initialized;  // set to true after the first spike from the synapse passes all the parameters

    double ltp_trace;
    double ltd_trace;

    std::deque< histentry_sh > ltp_ltd_history_;
//    std::deque< histentry_sh > weight_history_;

    //! Synapse specific parameters, so they can't be global members...
    double tau_ltp_;
    double tau_ltd_;

    double Tp_max_;
    double Td_max_;

    double eta_ltp_;  //!  Activation rate LTP trace
    double eta_ltd_;  //!  Activation rate LTD trace


    //! This is a local copy of the synaptic weight, mean for easier recording and debugging.
    //! Initialized to nan during construction, it's reinitialized/synced during the first query from
    //! the synapse with the real initial synaptic weight. Following this, the local weight is updated
    //! during every reward window locally so that we can keep track of the weight's evolution even if
    //! there are no spikes for a longer period.
    double local_weight_copy;
    double local_learn_rate_copy;

//    connectionDataStruct(double tau_ltp, double tau_ltd, double Tp_max, double Td_max, double eta_ltp, double eta_ltd):
    connectionDataStruct():
      initialized( false ),
      ltp_trace ( 0.0 ),
      ltd_trace ( 0.0 ),
      tau_ltp_( 2000. ),
      tau_ltd_( 1000. ),
      Tp_max_( 0.95 ),
      Td_max_( 1. ),
      eta_ltp_( 1. ),
      eta_ltd_( 1. ),
      local_weight_copy( inf_nan ),
      local_learn_rate_copy( 0.0 )
    {
    }
  };

  typedef std::map < index, connectionDataStruct* > connectionDataMap;


  void write_shouval_history( Time const& t_sp, index source_gid, double rate_pre, double rate_post, std::ostringstream &msg );
  void init_shouval_buffers();
  void get_status( DictionaryDatum& d ) const;
  void set_status( const DictionaryDatum& d );


  /**
 * For now, this function returns the trace for a single (first) synapse.
 * Later on it can be modified to compute the average over all.
 * @return
 */
  inline double
  get_ltp_trace_single() const
  {
    for (unsigned gid = 1; gid < 1000; ++gid)
    {
      if (connMap.find(gid) != connMap.end())
      {
        return connMap.at(gid)->ltp_trace;
      }
    }
    return 0.0;
  }

  inline double
  get_ltd_trace_single() const
  {
    for (unsigned gid = 1; gid < 1000; ++gid)
    {
      if (connMap.find(gid) != connMap.end())
      {
        return connMap.at(gid)->ltd_trace;
      }
    }
    return 0.0;
  }

  // TODO this can be made more efficient by keeping track of a running mean
  inline double
  get_mean_weight() const
  {
    int cnt = 0;
    double mean = 0;
    for ( auto it : connMap )
    {
       if ( it.second->local_weight_copy != inf_nan )
       {
         mean += it.second->local_weight_copy;
         cnt += 1;
       }
    }

    if ( cnt )
    {
      return mean / cnt;
    }
    return 0.0;
  }

  connectionDataMap get_connection_map()
  {
    return connMap;
  }


private:
//  std::map< index, std::vector< histentry_sh > > ltd_history_;
//  std::map< index, std::vector< histentry_sh > > ltp_history_;

//  unsigned int ltp_current;
  static const double inf_nan;
  unsigned int buffer_len;

  connectionDataMap connMap;

//  double tau_ltp_;
//  double tau_ltd_;
//
//  double Tp_max_;
//  double Td_max_;
//
//  double eta_ltp_;  //!  Activation rate LTP trace
//  double eta_ltd_;  //!  Activation rate LTD trace

  double T_tr_;  //!  Duration of refractory period for traces following neuromodulator presentation
  double T_reward_;  //!  duration of reward window

  std::vector< double > reward_times;

  /**
   *
   * @param t_ms
   * @return  0 - normal, 1 - reward, 2 - refractory
   */
  int
  get_window_type( const double t_ms )
  {
    for ( auto it_rew_time : reward_times )
    {
      if ( ( it_rew_time <= t_ms ) and ( t_ms < it_rew_time + T_reward_ ) )
      {
        return 1;
      }
      if ( ( it_rew_time + T_reward_ <= t_ms ) and ( t_ms < it_rew_time + T_reward_ + T_tr_ ) )
      {
        return 2;
      }
    }
    return 0;
  }

};


} // of namespace
#endif
