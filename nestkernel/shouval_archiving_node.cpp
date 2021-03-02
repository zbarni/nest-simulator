/*
 *  shouval_archiving_node.cpp
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

#include "shouval_archiving_node.h"

// Includes from nestkernel:
#include "kernel_manager.h"

// Includes from sli:
#include "dictutils.h"

namespace nest
{
const double nest::Shouval_Archiving_Node::inf_nan = std::numeric_limits<double>::max();


nest::Shouval_Archiving_Node::Shouval_Archiving_Node()
  : Archiving_Node()
  , buffer_len( 0 )
  , T_tr_( 25.0 )
  , T_reward_( 25.0 )
{
//  std::cout << "Constr Sh ArchivingNode " << this->get_gid() << std::endl << std::flush;
}

nest::Shouval_Archiving_Node::Shouval_Archiving_Node( const Shouval_Archiving_Node& n )
  : Archiving_Node( n )
  , buffer_len( n.buffer_len )
  , T_tr_( n.T_tr_ )
  , T_reward_( n.T_reward_ )
  , reward_times( n.reward_times )
{
//  std::cout << "CopyConstr Sh ArchivingNode " << this->get_gid() << std::endl << std::flush;
}

nest::Shouval_Archiving_Node::~Shouval_Archiving_Node()
{
  for ( auto it : connMap )
  {
    delete it.second;
  }
  connMap.clear();
}


void
nest::Shouval_Archiving_Node::init_shouval_buffers()
{
  buffer_len = kernel().connection_manager.get_max_delay() + 1;
}

void
nest::Shouval_Archiving_Node::get_status( DictionaryDatum& d ) const
{
  Archiving_Node::get_status( d );

  def<double>(d, names::T_tr, T_tr_);
  def<double>(d, names::T_reward, T_reward_);
  def<std::vector<double>>(d, names::reward_times, reward_times);
}

void
nest::Shouval_Archiving_Node::set_status( const DictionaryDatum& d )
{
  Archiving_Node::set_status( d );

  // We need to preserve values in case invalid values are set
  double new_T_tr = T_tr_;
  double new_T_reward = T_reward_;
  std::vector< double > new_reward_times = reward_times;

  updateValue<double>(d, names::T_tr, new_T_tr );
  updateValue<double>(d, names::T_reward, new_T_reward );
  updateValue< std::vector< double > >(d, names::reward_times, new_reward_times );

//  std::cout << "SET_STATUS " << this->get_gid() << " before " << T_tr_ << " after " << new_T_tr << std::endl << std::flush;
  T_tr_ = new_T_tr;
  T_reward_ = new_T_reward;
  reward_times = new_reward_times;
}

double
nest::Shouval_Archiving_Node::shouval_update_weight( double t_start, double t_stop,
                                                     double tau_ltp, double tau_ltd,
                                                     double Tp_max, double Td_max,
                                                     double eta_ltp, double eta_ltd,
                                                     double weight, double learn_rate,
                                                     index sender_gid )
{
//  std::cout << "[update weight] called for ["<<t_start<<","<<t_stop<<"] with w="<<weight << std::endl;

  if ( connMap.find( sender_gid ) == connMap.end())
  {
    return weight;
  }

  connectionDataStruct *conn_ptr = connMap[ sender_gid ];

  // initialize local parameters if hasn't been done before
  if ( !conn_ptr->initialized )
  {
    conn_ptr->initialized = true;
    conn_ptr->tau_ltp_ = tau_ltp;
    conn_ptr->tau_ltd_ = tau_ltd;
    conn_ptr->Tp_max_ = Tp_max;
    conn_ptr->Td_max_ = Td_max;
    conn_ptr->eta_ltp_ = eta_ltp;
    conn_ptr->eta_ltd_ = eta_ltd;
    conn_ptr->local_weight_copy = weight;
    conn_ptr->local_learn_rate_copy = learn_rate;
  }

  // if no history is stored, return weight
  if ( conn_ptr->ltp_ltd_history_.empty() )
  {
    return weight;
  }

  const double dt = nest::Time::get_resolution().get_ms();

  for ( histentry_sh h = conn_ptr->ltp_ltd_history_.front(); h.t_ <= t_stop and !conn_ptr->ltp_ltd_history_.empty();
        h = conn_ptr->ltp_ltd_history_.front() )
  {
//    std::cout << "[update weight] processing h_: " << h.t_ << std::endl;
    // remove item from front, will not be needed anymore
    conn_ptr->ltp_ltd_history_.pop_front();

    if ( t_start <= h.t_ )
    {
      // if local weight has been initialized before, it has been updated since and has to correct value at this time;
      // we can therefore avoid computing it again and just used the stored value
      if ( h.w_ != inf_nan )
      {
        weight = h.w_;
//        std::cout << "[update weight] valid h @: " << h.t_<<" ms, using precomputed weight: "<< weight << std::endl;
      }
      else  // update weight
      {
        double del_w = learn_rate * ( h.ltp_trace_ - h.ltd_trace_ ) * dt;
        weight += del_w;
//        std::cout << "[update weight] valid h @: " << h.t_<<" ms, new weight: "<< weight << std::endl;
      }
    }
  }
//  std::cout << "[update weight] returning new weight: " << weight << std::endl;
  return weight;
}

/**
 *
 * @param t_sp
 * @param source_gid
 * @param rate_pre
 * @param rate_post
 * @param msg
 */
void
nest::Shouval_Archiving_Node::write_shouval_history( Time const& t_sp, index source_gid,
                                                     double rate_pre, double rate_post, std::ostringstream &msg )
{
  const double t_ms = t_sp.get_ms();
  const double dt = nest::Time::get_resolution().get_ms();

  connectionDataStruct *conn_ptr;

  //! resize / initalize array for gid_pre when encountered for the first time
  if ( connMap.find( source_gid ) == connMap.end() )
  {
    conn_ptr = new connectionDataStruct();
    connMap[ source_gid ] = conn_ptr;

    msg << "\t\t|\tCreating connectionDataStruct() object for GID (pre) " << source_gid << std::endl;
  }
  else
  {
    conn_ptr = connMap[ source_gid ];
  }

  //! only update values after the first spikes, otherwise everything is still 0
  if ( conn_ptr->initialized )
  {
    //! hebbian term
    double h = rate_pre * rate_post;

    //! update traces
    double tmp_cur_ltp = conn_ptr->ltp_trace;
    double tmp_cur_ltd = conn_ptr->ltd_trace;

    // 0 - normal, 1 - reward, 2 - refractory
    int current_window = get_window_type( t_ms );

    if ( current_window < 2 )
    {
      double del_trace_ltp = ( -tmp_cur_ltp + conn_ptr->eta_ltp_ * h * ( conn_ptr->Tp_max_ - tmp_cur_ltp ) );
      double del_trace_ltd = ( -tmp_cur_ltd + conn_ptr->eta_ltd_ * h * ( conn_ptr->Td_max_ - tmp_cur_ltd ) );
      conn_ptr->ltp_trace += del_trace_ltp / conn_ptr->tau_ltp_ * dt;
      conn_ptr->ltd_trace += del_trace_ltd / conn_ptr->tau_ltd_ * dt;

      // store traces in reward window and update local weights if initialized
      if ( current_window == 1 )
      {
        msg << "\n\t\t|\t ... IN REWARD WINDOW \n";

        // update local weight copy
        if ( conn_ptr->local_weight_copy != inf_nan )
        {
          double del_w = conn_ptr->local_learn_rate_copy * ( conn_ptr->ltp_trace - conn_ptr->ltd_trace ) * dt;
          conn_ptr->local_weight_copy += del_w;
        }

        // store trace history entry
        conn_ptr->ltp_ltd_history_.push_back(
          histentry_sh( t_ms, conn_ptr->ltp_trace, conn_ptr->ltd_trace, conn_ptr->local_weight_copy, 0 ) );

        msg << "\n\t\t|\t ... inserted new history entry with t_ms = " << t_ms << " ... \n";
      }
    }
    else
    {
      conn_ptr->ltp_trace = 0.0;
      conn_ptr->ltd_trace = 0.0;
    }

    msg << "\t\t|\tUpdating trace for GID (pre) " << source_gid << " @t: " << t_ms << " ms"
        << "\n\t\t|\t current trace value " << tmp_cur_ltp << " ->  " << conn_ptr->ltp_trace
        //      << "\n\t current weight value " << tmp_weight << " ->  " << conn_ptr->weight
        << "\n\t\t|\t window type " << current_window << " len(reward_times): " << reward_times.size()
        << "\n\t\t|\t rate pre " << rate_pre << " and rate post " << rate_post << "\t hebbian " << h << std::endl
        << std::flush;
  }
  else
  {
    msg << "\t\t|\tNot yet updating trace for GID (pre) " << source_gid << " @t " << t_ms<< " ms; hasn't spiked yet!\n";
  }
}


} // of namespace nest
