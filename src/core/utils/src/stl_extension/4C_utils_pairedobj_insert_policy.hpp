// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_UTILS_PAIREDOBJ_INSERT_POLICY_HPP
#define FOUR_C_UTILS_PAIREDOBJ_INSERT_POLICY_HPP

#include "4C_config.hpp"

#include "4C_utils_exceptions.hpp"

#include <algorithm>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <ostream>
#include <sstream>
#include <utility>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Core::Gen
{
  /// binary function
  template <typename T>
  inline bool pair_comp(const T& i, const T& j)
  {
    return i.first < j.first;
  }

  /*--------------------------------------------------------------------------*/
  template <typename Key, typename T>
  class DefaultInsertPolicy
  {
    typedef DefaultInsertPolicy<Key, T> class_type;

   protected:
    typedef std::pair<Key, T> pair_type;
    typedef std::vector<pair_type> pairedvector_type;
    typedef typename pairedvector_type::iterator iterator;
    typedef typename pairedvector_type::const_iterator const_iterator;

   public:
    iterator begin(const iterator& first) { return first; }

    const_iterator begin(const const_iterator& first) const { return first; }

    iterator end(const iterator& last) { return last; }

    const_iterator end(const const_iterator& last) const { return last; }

    iterator find(const Key k, pairedvector_type& data, size_t entries)
    {
      iterator last = data.begin() + entries;
      for (iterator it = data.begin(); it != last; ++it)
      {
        if (it->first == k) return it;
      }
      return last;
    }

    const_iterator find(const Key k, const pairedvector_type& data, size_t entries) const
    {
      const_iterator last = data.begin() + entries;
      for (const_iterator it = data.begin(); it != last; ++it)
      {
        if (it->first == k) return it;
      }
      return last;
    }

    T& get(const Key k, pairedvector_type& data, size_t& entries)
    {
      iterator last = data.begin() + entries;
      iterator it = find(k, data, entries);
      if (it != last) return it->second;

      if (entries >= data.size()) throw std::length_error("Pairedvector::operator[]");

      ++entries;
      last->first = k;
      return last->second;
    }

    T& at(const Key k, pairedvector_type& data, size_t entries)
    {
      iterator last = data.begin() + entries;
      iterator it = find(k, data, entries);
      if (it == last) FOUR_C_THROW("Pairedvector::at(): invalid key");

      return it->second;
    }

    const T& at(const Key k, const pairedvector_type& data, size_t entries) const
    {
      const_iterator last = data.begin() + entries;
      const_iterator it = find(k, data, entries);
      if (it == last) FOUR_C_THROW("Pairedvector::at(): invalid key");

      return it->second;
    }

    T& operator()(const Key k, pairedvector_type& data, size_t& entries)
    {
      return get(k, data, entries);
    }

    T& repetitive_access(const Key k, const int rep_count, pairedvector_type& data, size_t& entries)
    {
      return get(k, data, entries);
    }

    size_t complete(pairedvector_type& unique_vec, const size_t unique_entries)
    {
      return unique_entries;
    }

    void swap(class_type& x) { /* empty */ }

    void clone(const class_type& source) { /* empty */ }

    /** @brief return a constant number of offset values
     *
     *   The default strategy does not use this feature. See quick_insert_policy
     *   for an example.
     *
     */
    static size_t capacity_offset() { return 0; }
  };

  /*--------------------------------------------------------------------------*/
  template <typename Key, typename T>
  class QuickInsertPolicy : public DefaultInsertPolicy<Key, T>
  {
    typedef QuickInsertPolicy<Key, T> class_type;

   protected:
    typedef DefaultInsertPolicy<Key, T> base_type;
    typedef typename base_type::pair_type pair_type;
    typedef typename base_type::pairedvector_type pairedvector_type;
    typedef typename base_type::iterator iterator;
    typedef typename base_type::const_iterator const_iterator;

   public:
    /** @brief return a constant number of offset values
     *
     *  This number will be added to a new capacity and can be used for
     *  sentinel values for example.
     *
     */
    static size_t capacity_offset() { return 1; }

    T& get(const Key k, pairedvector_type& data, size_t& entries)
    {
      iterator last = data.begin() + entries;
      iterator it = find(k, data, entries);
      if (it != last) return it->second;

      if (entries + capacity_offset() >= data.size())
        throw std::length_error("Pairedvector::get()");

      ++entries;
      return last->second;
    }

    /// copy from the base class
    T& at(const Key k, pairedvector_type& data, size_t entries)
    {
      iterator last = data.begin() + entries;
      iterator it = find(k, data, entries);
      if (it == last) FOUR_C_THROW("Pairedvector::at(): invalid key");

      return it->second;
    }

    /// copy from the base class
    const T& at(const Key k, const pairedvector_type& data, size_t entries) const
    {
      const_iterator last = data.begin() + entries;
      const_iterator it = find(k, data, entries);
      if (it == last) FOUR_C_THROW("Pairedvector::at(): invalid key");

      return it->second;
    }

    /// copy from the base class
    T& operator()(const Key k, pairedvector_type& data, size_t& entries)
    {
      return get(k, data, entries);
    }

    /** @brief Quick linear search routine
     *
     *  The unnecessary bound check of the for loop is avoided by introduction
     *  of a sentinel value at the one past the last position.
     *
     *  @param[in] k       Unique key value.
     *  @param[in] data    Internally stored data.
     *  @param[in] entries Number of current entries.
     *  @return A iterator to the pair with the given key is returned. If the
     *          key is not a part of the stored data, an iterator to one past
     *          the last position will be returned.
     *
     */
    iterator find(const Key k, pairedvector_type& data, size_t entries)
    {
      // set sentinel values
      (data.begin() + entries)->first = k;

      iterator it = data.begin();
      for (;;)
      {
        if (it->first == k) return it;
        ++it;
      }

      exit(EXIT_FAILURE);
    }

    /** @brief Quick linear search routine (const version)
     *
     *  The unnecessary bound check of the for loop is avoided by introduction
     *  of a sentinel value at the one past the last position. Therefore a
     *  technically const_cast becomes necessary. Nevertheless, the real data
     *  set is not touched.
     *
     *  @param[in] k       Unique key value.
     *  @param[in] data    Internally stored data.
     *  @param[in] entries Number of current entries.
     *  @return A const_iterator to the pair with the given key is returned. If
     *          the key is not a part of the stored data, an iterator to one past
     *          the last position will be returned.
     *
     */
    const_iterator find(const Key k, const pairedvector_type& data, size_t entries) const
    {
      // set sentinel values
      pairedvector_type& mutable_data = const_cast<pairedvector_type&>(data);
      (mutable_data.begin() + entries)->first = k;

      const_iterator cit = data.begin();
      for (;;)
      {
        if (cit->first == k) return cit;
        ++cit;
      }

      exit(EXIT_FAILURE);
    }

    /** @brief repetitive access of the paired_vector entries in the same order
     *
     *  This method is very helpful, if a big paired_vector structure must be
     *  accessed repetitively in the same order. In this case only the first
     *  repetition (i.e. rep_count == 0) uses the internal linear search
     *  routine, while all following access repetitions (i.e. rep_count > 0)
     *  use the internal id-map under the assumption of an constant access
     *  pattern.
     *
     *  \note This method can not be used, if the access pattern changes during
     *  successive repetitions!
     *
     *  @param[in] k         Unique identification key of the stored data.
     *  @param[in] rep_count Repetition counter.
     *  @param[in] data      Internal stored data set.
     *  @param[in] entries   Updated number of entries.
     *  @return The value corresponding to the given key.
     *
     */
    T& repetitive_access(const Key k, const int rep_count, pairedvector_type& data, size_t& entries)
    {
#ifdef FOUR_C_ENABLE_ASSERTIONS
      FOUR_C_ASSERT(rep_count < -1,
          "The repetition counter is not allowed to "
          "be smaller than -1!");
#endif

      switch (rep_count)
      {
        case -1:
        {
          return get(k, data, entries);
        }
        case 0:
        {
#ifdef FOUR_C_ENABLE_ASSERTIONS
          FOUR_C_ASSERT(prev_repcount_ != -1 and prev_repcount_ != rep_count and !isfilled_,
              "At the end of each repetition Complete "
              "must be called!");
#endif
          init_id_map();
          isfilled_ = false;
          return first_access(k, data, entries);
        }
        default:
        {
#ifdef FOUR_C_ENABLE_ASSERTIONS
          FOUR_C_ASSERT(prev_repcount_ != rep_count and !isfilled_,
              "At the end of each repetition Complete "
              "must be called!");
#endif
          isfilled_ = false;
          prev_repcount_ = rep_count;
          return subsequent_access(k, data, entries);
        }
      }
    }

    size_t complete(pairedvector_type& unique_vec, const size_t unique_entries)
    {
      if (isfilled_) return unique_entries;

      // These things must be done only once after the first complete call,
      // i.e. repetition counter equal to zero.
      if (prev_repcount_ == 0)
      {
        id_map_size_ = id_map_.size();
        if (id_index_ != id_map_size_)
        {
          FOUR_C_THROW("Size (== {}) <--> count (== {}) mismatch!\n", id_map_.size(), id_index_);
        }
      }

      isfilled_ = true;
      // reset id-index after each single complete call
      id_index_ = 0;

      return unique_entries;
    }

    void clone(const class_type& source)
    {
      base_type::clone(source);

      id_map_ = source.id_map_;
      id_map_size_ = source.id_map_size_;
      id_index_ = source.id_index_;

      prev_repcount_ = source.prev_repcount_;
      isfilled_ = source.isfilled_;
    }

   private:
    /** @brief Initiate the internally stored ID map
     *
     *  The old map is cleared, while the old size is used as an estimate for
     *  the necessary storage space.
     *
     */
    void init_id_map()
    {
      // nothing to do, if it is a successive access in the first repetition
      // Note: The _isfilled check guarantees that the maps are cleared, even if
      // only one GP had been considered in the previous attempt.
      if (prev_repcount_ == 0 and not isfilled_) return;

      id_map_.clear();
      id_map_.reserve(id_map_size_);
      prev_repcount_ = 0;
    }

    T& first_access(const Key k, pairedvector_type& data, size_t& entries)
    {
      iterator last = data.begin() + entries;
      iterator it = find(k, data, entries);

      // By using push_back we follow the std::vector policy if the set
      // capacity is exceeded. See previous reserve call in init_id_map().
      id_map_.push_back(it - data.begin());
      ++id_index_;

      if (it == last)
      {
        if (entries + capacity_offset() >= data.size())
          throw std::length_error("quick_insert_policy::get()");

        ++entries;
      }

      return it->second;
    }

    T& subsequent_access(const Key k, pairedvector_type& data, size_t& entries)
    {
      pair_type& curr_pair = data[get_id()];
      if (curr_pair.first != k)
      {
        FOUR_C_THROW(
            "quick_insert_policy::subsequent_access(): Mismatch of stored "
            "key and given key: {} <--> {}",
            curr_pair.first, k);
      }

      return curr_pair.second;
    }

    size_t get_id() { return id_map_[id_index_++]; }

   private:
    std::vector<size_t> id_map_;
    size_t id_map_size_ = 512;
    size_t id_index_ = 0;

    int prev_repcount_ = -1;
    bool isfilled_ = false;
  };

  /*--------------------------------------------------------------------------*/
  template <typename Key, typename T>
  class InsertAndSortPolicy : public DefaultInsertPolicy<Key, T>
  {
    typedef InsertAndSortPolicy<Key, T> class_type;

   protected:
    typedef DefaultInsertPolicy<Key, T> base_type;
    typedef typename base_type::pair_type pair_type;
    typedef typename base_type::pairedvector_type pairedvector_type;
    typedef typename base_type::iterator iterator;
    typedef typename base_type::const_iterator const_iterator;

   public:
    InsertAndSortPolicy()
        : non_unique_vec_(0),
          non_unique_entries_(0),
          dyn_max_allowed_capacity_(CONST_MAX_ALLOWED_CAPACITY),
          isfilled_(true),
          max_value_(0.0) {};

    iterator begin(const iterator& first)
    {
      throw_if_not_filled(__LINE__, __FUNCTION__);
      return base_type::begin(first);
    }

    const_iterator begin(const const_iterator& first) const
    {
      throw_if_not_filled(__LINE__, __FUNCTION__);
      return base_type::begin(first);
    }

    iterator end(const iterator& last)
    {
      throw_if_not_filled(__LINE__, __FUNCTION__);
      return base_type::end(last);
    }

    const_iterator end(const const_iterator& last) const
    {
      throw_if_not_filled(__LINE__, __FUNCTION__);
      return base_type::end(last);
    }

    iterator find(const Key k, pairedvector_type& data, size_t entries)
    {
      throw_if_not_filled(__LINE__, __FUNCTION__);
      return base_type::find(k, data, entries);
    }

    const_iterator find(const Key k, const pairedvector_type& data, size_t entries) const
    {
      throw_if_not_filled(__LINE__, __FUNCTION__);
      return base_type::find(k, data, entries);
    }

    T& get(const Key k, pairedvector_type& data, size_t& entries)
    {
      throw_if_not_filled(__LINE__, __FUNCTION__);
      return base_type::get(k, data, entries);
    }

    T& at(const Key k, pairedvector_type& data, size_t entries)
    {
      throw_if_not_filled(__LINE__, __FUNCTION__);
      return base_type::at(k, data, entries);
    }

    const T& at(const Key k, const pairedvector_type& data, size_t entries) const
    {
      throw_if_not_filled(__LINE__, __FUNCTION__);
      return base_type::at(k, data, entries);
    }

    T& repetitive_access(const Key k, const int rep_count, pairedvector_type& data, size_t& entries)
    {
      return this->operator()(k, data, entries);
    }

    /** @brief Insert values into the non-unique vector
     *
     *  @note This routine makes it necessary that you call complete before you can
     *  access the inserted values.
     *  @note Values with an absolute value smaller than the relative machine
     *  precision are ignored.
     *
     *  @param[in] k       Unique key value.
     *  @param[in] data    Unique paired data vector (unused at this point)
     *  @param[in] entries Number of unique data in the data vector.
     *  @return Reference to a new slot to insert data.
     *
     */
    T& operator()(const Key k, pairedvector_type& data, size_t& entries)
    {
      if (non_unique_entries_ == dyn_max_allowed_capacity_ and
          (not set_dynamic_max_allowed_capacity(2 * entries)))
      {
        entries = mid_complete(data, entries);
      }

      if (isfilled_)
      {
        isfilled_ = false;
        init_capacity();
        ++non_unique_entries_;
      }

      pair_type* pair_ptr = non_unique_vec_.data() + current_id();

      set_max_value(pair_ptr->second);

      /* If reasonable values have been inserted, we check the capacity and go
       * to the next entry. Otherwise, we stick to the current entry. */
      if (std::abs(pair_ptr->second) > get_relative_machine_precision())
      {
        increase_capacity();
        ++non_unique_entries_;
        pair_ptr = non_unique_vec_.data() + current_id();
      }

      pair_ptr->first = k;
      return pair_ptr->second;
    }

    inline bool is_filled() const { return isfilled_; }

    size_t complete(pairedvector_type& unique_vec, const size_t unique_entries)
    {
      if (isfilled_) return unique_entries;

      const size_t new_unique_entries = group_and_merge(unique_vec, unique_entries);

      final_post_complete();

      return new_unique_entries;
    }

    void swap(class_type& x)
    {
      base_type::swap(x);

      const bool my_isfilled = isfilled_;
      isfilled_ = x.isfilled_;
      x.isfilled_ = my_isfilled;

      const size_t my_non_unique_entries = non_unique_entries_;
      non_unique_entries_ = x.non_unique_entries_;
      x.non_unique_entries_ = my_non_unique_entries;

      const size_t my_dyn_max_allowed_capacity = dyn_max_allowed_capacity_;
      dyn_max_allowed_capacity_ = x.dyn_max_allowed_capacity_;
      x.dyn_max_allowed_capacity_ = my_dyn_max_allowed_capacity;

      const double my_max_value = max_value_;
      max_value_ = x.max_value_;
      x.max_value_ = my_max_value;

      non_unique_vec_.swap(x.non_unique_vec_);
    }

    void clone(const class_type& source) { throw_if_not_filled(__LINE__, __FUNCTION__); }

   private:
    /** @brief Internal complete call
     *
     *  This is a complete call done during the a standard insert call into the
     *  non-unique vector. Its purpose is to group and merge the non-unique data
     *  as soon as some bound for the size of the non-unique vector is reached.
     *  In contrast to a default complete call the already allocated memory of
     *  the non-unique vector is not freed afterwards, since the insertion of
     *  new values is not yet finished.
     *
     *  @param[in/out] unique_vec     This object will contain the new
     *                                completed unique data.
     *  @param[in]     unique_entries Previous/old number of unique data.
     *  @return The number of the unique data in the end of this routine.
     *
     */
    size_t mid_complete(pairedvector_type& unique_vec, const size_t unique_entries)
    {
      if (isfilled_) return unique_entries;

      const size_t new_unique_entries = group_and_merge(unique_vec, unique_entries);

      mid_post_complete();

      return new_unique_entries;
    }

    /** @brief Throw a run time error if any method is called which
     *  presumes a filled state
     *
     *  @param[in] linenumber  Line number in the calling function.
     *  @param[in] functname   Function name of the calling function.
     *
     */
    inline void throw_if_not_filled(int linenumber, const std::string& functname) const
    {
      if (not isfilled_)
      {
        FOUR_C_THROW("LINE {} in {}: access denied! call complete first.", linenumber, functname);
      }
    }

    /** @brief Group all entries in the %_non_unique_vec and insert the
     *  merged values of each group into the %unique_vec
     *
     *  @param[out] unique_vec     Vector containing the current unique entries.
     *                             This vector will be modified by inserting
     *                             the entries collected in the non-unique vector.
     *  @param[in]  unique_entries Number of old unique entries.
     *  @return Number of the new unique entries contained in unique_vec.
     *
     */
    size_t group_and_merge(pairedvector_type& unique_vec, const size_t unique_entries)
    {
      // get total number of entries and merge all entries temporal in the
      // _non_unique_vec member
      const size_t merged_entries = unique_entries + non_unique_entries_;
      if (non_unique_vec_.size() < merged_entries + 1)
        non_unique_vec_.resize(merged_entries + 1, pair_type());

      std::copy(unique_vec.begin(), unique_vec.begin() + unique_entries,
          non_unique_vec_.begin() + non_unique_entries_);

      // group the entries
      const size_t num_grps = group_data(non_unique_vec_.begin(), merged_entries);

      if (unique_vec.size() < num_grps) unique_vec.resize(num_grps);

      // copy and merge all non-unique entries into one unique vector
      const typename pairedvector_type::iterator last_result = merge_group_data(
          non_unique_vec_.begin(), non_unique_vec_.begin() + merged_entries, unique_vec.begin());

      const size_t corrected_length = static_cast<size_t>(last_result - unique_vec.begin());

      // reset the remaining part of the unique vector
      if (corrected_length < unique_entries)
        std::fill(last_result, unique_vec.begin() + unique_entries, pair_type());

      return corrected_length;
    }

    /** @brief reset stuff after an internal complete call
     *
     *  In contrast to the %final_post_complete method we do not free the
     *  allocated memory for the %_non_unique_vec, instead we just clear the
     *  content. The %_max_value and the current necessary dynamic
     *  maximal allowed capacity are also kept.
     *
     */
    void mid_post_complete()
    {
      std::fill(
          non_unique_vec_.begin(), non_unique_vec_.begin() + non_unique_entries_, pair_type());

      if (_max_capacity < non_unique_entries_) _max_capacity = non_unique_entries_;

      non_unique_entries_ = 0;
      isfilled_ = true;
    }

    /** @brief reset all class members at the end of a final externally
     *  initialized complete call
     *
     *  The allocated memory for the %_non_unique_vec is freed. Furthermore,
     *  all class members which are only important during one call are reset.
     *
     */
    void final_post_complete()
    {
      non_unique_vec_.clear();

      if (_max_capacity < non_unique_entries_) _max_capacity = non_unique_entries_;

      dyn_max_allowed_capacity_ = CONST_MAX_ALLOWED_CAPACITY;
      max_value_ = 0.0;
      non_unique_entries_ = 0;
      isfilled_ = true;
    }

    /** @brief Group data into clusters with same KEY values
     *
     *  Group data in the given interval by creating clusters of pairs with the
     *  same KEY value.
     *
     *  @param[in] sbegin  iterator pointing at the first element of the interval
     *  @return Number of groups.
     *
     */
    size_t group_data(typename pairedvector_type::iterator sbegin, const size_t num_entries) const
    {
      switch (num_entries)
      {
        case 0:
          return 0;
        case 1:
          return 1;
        case 2:
        {
          if (sbegin->first == (sbegin + 1)->first) return 1;
          return 2;
        }
        default:
          return group_big_data(sbegin, sbegin + num_entries);
      }
    }

    /** @brief Group data into clusters with same KEY values
     *
     *  Group data in the given interval by creating clusters of pairs with the
     *  same KEY value. No extra memory is allocated since the grouping works
     *  by a simple swap routine.
     *
     *  @param[in] sbegin  iterator pointing at the first element of the interval
     *  @param[in] slast   iterator pointing one past the last element of the interval
     *  @return Number of groups.
     *
     */
    size_t group_big_data(typename pairedvector_type::iterator sbegin,
        const typename pairedvector_type::iterator slast) const
    {
      typename pairedvector_type::iterator result = sbegin;
      const typename pairedvector_type::iterator prev_last = slast - 1;

      size_t num_grps = 1;

      for (;;)
      {
        // set sentinel accordingly
        *slast = *result;
        sbegin = result + 1;

        for (;;)
        {
          pair_type& val = *sbegin;
          if (val.first == result->first)
          {
            if (sbegin == slast) break;

            // swap values
            (++result)->swap(val);
          }
          ++sbegin;
        }

        if (result == prev_last) break;

        ++result;
        ++num_grps;
      }

      return num_grps;
    }

    /** @brief Merge the values in each group
     *
     *  @pre The input values in between %sbegin and %slast must be grouped.
     *
     *  @param[in]  sbegin  constant iterator pointing to the first entry of
     *                      already grouped interval
     *  @param[in]  slast   iterator pointing to one past the last
     *                      entry of the already grouped interval (used as
     *                      sentinel value)
     *  @param[out] result  iterator pointing to the merged result interval
     *  @return Return the iterator pointing to one past the last entry of the
     *          result interval
     *
     */
    typename pairedvector_type::iterator merge_group_data(
        typename pairedvector_type::const_iterator sbegin,
        typename pairedvector_type::iterator slast,
        typename pairedvector_type::iterator result) const
    {
      if (sbegin == slast) return result;

      // add sentinel
      slast->first = ((slast - 1)->first + 1);

      *result = *sbegin;
      for (;;)
      {
        const std::pair<Key, T>& val = *(++sbegin);
        if (val.first == result->first)
          result->second += val.second;
        else
        {
          if (sbegin == slast) break;

          // erase/overwrite zero entries
          if (std::abs(result->second) > get_relative_machine_precision()) ++result;

          *result = val;
        }
      }

      // Kick the last element out if the absolute value of the last result is
      // less than the relative machine precision.
      return (std::abs(result->second) > get_relative_machine_precision() ? ++result : result);
    }

    /// Return the local id of the last element
    inline size_t current_id() const { return (non_unique_entries_ - 1); }

    /** @brief Set a new maximal allowed capacity of the non-unique vector
     *
     *  The maximal allowed capacity is only modified if the new dynamic bound
     *  value is larger as the current one.
     *
     *  @param[in] dyn_bound New dynamic bound value.
     *  @return If the threshold is modified, TRUE is returned, otherwise FALSE.
     *
     */
    bool set_dynamic_max_allowed_capacity(const size_t dyn_bound)
    {
      if (dyn_bound > dyn_max_allowed_capacity_)
      {
        dyn_max_allowed_capacity_ = dyn_bound;
        return true;
      }
      return false;
    }

    /** @brief Set maximal absolute value occurring during the current insertion
     *         interval, i.e. between the last and next complete call.
     *
     *  @param[in] curr_val This is the current value which is going to be inserted.
     *
     */
    void set_max_value(const double curr_val)
    {
      const double abs_curr_val = std::abs(curr_val);
      if (abs_curr_val > max_value_) max_value_ = abs_curr_val;
    }

    /** @brief Get the relative machine precision
     *
     *  Relative to the max value occurring during the current insertion
     *  interval. See set_max_value for more information.
     *
     */
    inline double get_relative_machine_precision() const { return max_value_ * MACHINE_PRECISION; }

    inline void init_capacity()
    {
      non_unique_vec_.resize(std::max(1, static_cast<int>(_max_capacity)), pair_type());
    }

    void increase_capacity()
    {
      const size_t curr_length = non_unique_vec_.size();
      if (non_unique_entries_ < curr_length) return;

      non_unique_vec_.resize(2 * curr_length, pair_type());
    }

    void print(const pairedvector_type& pvec, std::ostream& os = std::cout) const
    {
      os << "[Size = " << pvec.size() << "]\n"
         << "{ ";
      size_t c = 0;
      typename pairedvector_type::const_iterator first = pvec.begin();
      while (first != pvec.end())
      {
        const pair_type& pair = *first;
        os << "[" << c++ << "](" << pair.first << ", " << pair.second << ")";
        if (++first != pvec.end()) os << "; ";
      }
      os << " };\n" << std::flush;
    }

   private:
    /** temporal allocated vector containing the unsorted and therefore
     *  non-unique entries */
    pairedvector_type non_unique_vec_;

    /** @brief number of the inserted non-unique entries
     *
     *  Note that this number is in general not equal to the size of the
     *  _non_unique_vec, since memory is pre-allocated. */
    size_t non_unique_entries_;

    /** maximal necessary capacity for all paired-vector objects of the same
     *  types */
    static size_t _max_capacity;

    /// absolute machine precision
    static constexpr double MACHINE_PRECISION = 10 * std::numeric_limits<double>::epsilon();

    /// constant maximal allowed capacity
    static constexpr size_t CONST_MAX_ALLOWED_CAPACITY = 512;

    /** @brief Dynamically corrected maximal allowed capacity.
     *
     *  If the unique vector exceeds the constant capacity, the upper capacity
     *  bound must be increased to stay meaningful.
     *
     *  See the method %set_dynamic_max_allowed_capacity() for more information. */
    size_t dyn_max_allowed_capacity_;

    bool isfilled_;
    double max_value_;
  };
}  // namespace Core::Gen

template <typename Key, typename T>
size_t Core::Gen::InsertAndSortPolicy<Key, T>::_max_capacity = 0;



FOUR_C_NAMESPACE_CLOSE

#endif
