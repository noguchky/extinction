#ifndef Tron_Linq_hh
#define Tron_Linq_hh

#include <algorithm>
#include <numeric>
#include <map>
#include <set>
#include <string>
#include <sstream>
#include <vector>

namespace Tron {

  template <typename Key_t, typename Value_t>
  using Lookup_t = std::map<Key_t, std::vector<Value_t>>;

  namespace Linq {

    template <typename T>
    class Enumerable_t {
    public:
      using Iterable_t = T;
      using Iterator_t = typename Iterable_t::const_iterator;
      using Value_t    = typename Iterable_t::value_type;

    private:
      Iterable_t fIterable;

    public:
      Enumerable_t(const Enumerable_t& enumerable)
        : fIterable(enumerable.fIterable) {
      }

      Enumerable_t(const Iterable_t& iterable)
        : fIterable(iterable) {
      }

      Enumerable_t(const Value_t* begin, const Value_t* end)
        : fIterable(begin, end) {
      }

      template <typename Iterator_t>
      Enumerable_t(const Iterator_t begin, const Iterator_t end)
        : fIterable(begin, end) {
      }

      auto begin() const -> Iterator_t {
        return fIterable.begin();
      }

      auto end() const -> Iterator_t {
        return fIterable.end();
      }

      auto operator->() const -> Iterable_t* {
        return &fIterable;
      }

      auto IsEmpty() const -> bool {
        return begin() == end();
      }

      template <typename Predicate_t>
      auto Any(Predicate_t predicate) const -> bool {
        return std::any_of(begin(), end(), predicate);
      }

      template <typename Predicate_t>
      auto All(Predicate_t predicate) const -> bool {
        return std::all_of(begin(), end(), predicate);
      }

      auto Contains(const Value_t& value) const -> bool {
        return std::find(begin(), end(), value) != end();
      }

      template <typename Predicate_t,
                typename NewIterable_t = std::vector<Value_t>>
      auto Where(Predicate_t predicate) const -> Enumerable_t<NewIterable_t> {
        NewIterable_t values;
        std::copy_if(begin(), end(), std::back_inserter(values), predicate);
        return values;
      }

      template <typename Predicate_t,
                typename NewValue_t    = decltype((*((Predicate_t*)nullptr))(Value_t())),
                typename NewIterable_t = std::vector<NewValue_t>>
      auto Select(Predicate_t predicate) const -> Enumerable_t<NewIterable_t> {
        NewIterable_t values;
        std::transform(begin(), end(), std::back_inserter(values), predicate);
        return values;
      }

      template <typename Predicate_t>
      auto ForEach(Predicate_t predicate) const -> Enumerable_t {
        std::for_each(begin(), end(), predicate);
        return *this;
      }

      auto Sum() const -> Value_t {
        return std::accumulate(begin(), end(), Value_t());
      }

      template <typename Predicate_t,
                typename NewValue_t    = decltype((*((Predicate_t*)nullptr))(Value_t()))>
      auto Sum(Predicate_t predicate) const -> NewValue_t {
        std::vector<NewValue_t> values;
        std::transform(begin(), end(), std::back_inserter(values), predicate);
        return std::accumulate(values.begin(), values.end(), NewValue_t());
      }

      template <typename NewIterable_t = std::vector<Value_t>>
      auto Skip(unsigned int n) const -> Enumerable_t<NewIterable_t> {
        if (n >= Count()) {
          return NewIterable_t();
        }
        return NewIterable_t(std::next(begin(), n), end());
      }

      template <typename NewIterable_t = std::vector<Value_t>>
      auto Take(unsigned int n) const -> Enumerable_t<NewIterable_t> {
        if (n >= Count()) {
          return NewIterable_t(begin(), end());
        }
        return NewIterable_t(begin(), std::next(begin(), n));
      }

      template <typename Predicate_t,
                typename NewIterable_t = std::vector<Value_t>>
      auto SkipWhile(Predicate_t predicate) const -> Enumerable_t<NewIterable_t> {
        unsigned int n = 0;
        for (auto&& value : *this) {
          if (predicate(value)) {
            n++;
          } else {
            break;
          }
        }
        return Skip<NewIterable_t>(n);
      }

      template <typename Predicate_t,
                typename NewIterable_t = std::vector<Value_t>>
      auto TakeWhile(Predicate_t predicate) const -> Enumerable_t<NewIterable_t> {
        unsigned int n = 0;
        for (auto&& value : *this) {
          if (predicate(value)) {
            n++;
          } else {
            break;
          }
        }
        return Take<NewIterable_t>(n);
      }

      auto Min() const -> Value_t {
        auto minitr = std::min_element(begin(), end());
        return minitr == end() ? Value_t() : *minitr;
      }

      auto Min(Value_t defValue) const -> Value_t {
        auto minitr = std::min_element(begin(), end());
        return minitr == end() ? defValue : *minitr;
      }

      auto Max() const -> Value_t {
        auto maxitr = std::max_element(begin(), end());
        return maxitr == end() ? Value_t() : *maxitr;
      }

      auto Max(Value_t defValue) const -> Value_t {
        auto maxitr = std::max_element(begin(), end());
        return maxitr == end() ? defValue : *maxitr;
      }

      auto First() const -> Value_t {
        return *begin();
      }

      auto FirstOrDefault() const -> Value_t {
        return IsEmpty() ? Value_t() : First();
      }

      auto Last() const -> Value_t {
        return *(end() - 1);
      }

      auto LastOrDefault() const -> Value_t {
        return IsEmpty() ? Value_t() : Last();
      }

      auto Count() const -> unsigned int {
        return std::distance(begin(), end());
      }

      template <typename Predicate_t>
      auto Count(Predicate_t predicate) const -> unsigned int {
        return std::count_if(begin(), end(), predicate);
      }

      auto Average() const -> Value_t {
        return IsEmpty() ? Sum() / Count() : Value_t();
      }

      auto Distinct() const -> Enumerable_t<std::set<Value_t>> {
        return std::set<Value_t>(begin(), end());
      }

      template <typename OValue_t,
                typename NewIterable_t = std::vector<std::pair<Value_t, OValue_t>>>
      auto Join(const OValue_t* begin) -> Enumerable_t<NewIterable_t> {
        NewIterable_t joined;
        auto it1 = this->begin();
        auto it2 = begin;
        for (; it1 != this->end(); ++it1, ++it2) {
          joined.push_back({ *it1, *it2 });
        }
        return joined;
      }

      template <typename OIterator_t,
                typename NewIterable_t = std::vector<std::pair<Value_t, typename OIterator_t::value_type>>>
      auto Join(const OIterator_t& begin) -> Enumerable_t<NewIterable_t> {
        NewIterable_t joined;
        auto it1 = this->begin();
        auto it2 = begin;
        for (; it1 != this->end(); ++it1, ++it2) {
          joined.push_back({ *it1, *it2 });
        }
        return joined;
      }

      template <typename V>
      auto To() const -> V {
        return V(begin(), end());
      }

      auto ToVector() const -> std::vector<Value_t> {
        return std::vector<Value_t>(begin(), end());
      }

      auto ToSet() const -> std::set<Value_t> {
        return std::set<Value_t>(begin(), end());
      }

      template <typename   KeyPredicate_t,
                typename ValuePredicate_t,
                typename NewKey_t   = decltype((*((  KeyPredicate_t*)nullptr))(Value_t())),
                typename NewValue_t = decltype((*((ValuePredicate_t*)nullptr))(Value_t()))>
      auto ToMap(KeyPredicate_t keyPredicate, ValuePredicate_t valuePredicate) const -> std::map<NewKey_t, NewValue_t> {
        std::map<NewKey_t, NewValue_t> map;
        for (auto&& value : *this) {
          NewKey_t   key      =   keyPredicate(value);
          NewValue_t newValue = valuePredicate(value);
          map[key] = newValue;
        }
        return map;
      }

      template <typename KeyPredicate_t,
                typename NewKey_t = decltype((*((KeyPredicate_t*)nullptr))(Value_t()))>
      auto ToMap(KeyPredicate_t keyPredicate) const -> std::map<NewKey_t, Value_t> {
        auto valuePredicate = [](Value_t value) -> Value_t { return value; };
        return ToMap<decltype(keyPredicate), decltype(valuePredicate), NewKey_t, Value_t>(keyPredicate, valuePredicate);
      }

      template <typename   KeyPredicate_t,
                typename ValuePredicate_t,
                typename NewKey_t   = decltype((*((  KeyPredicate_t*)nullptr))(Value_t())),
                typename NewValue_t = decltype((*((ValuePredicate_t*)nullptr))(Value_t()))>
      auto ToLookup(KeyPredicate_t keyPredicate, ValuePredicate_t valuePredicate) const -> Lookup_t<NewKey_t, NewValue_t> {
        Lookup_t<NewKey_t, NewValue_t> lookup;
        for (auto&& value : *this) {
          NewKey_t   key      =   keyPredicate(value);
          NewValue_t newValue = valuePredicate(value);
          lookup[key].push_back(newValue);
        }
        return lookup;
      }

      template <typename KeyPredicate_t,
                typename NewKey_t = decltype((*((KeyPredicate_t*)nullptr))(Value_t()))>
      auto ToLookup(KeyPredicate_t keyPredicate) const -> Lookup_t<NewKey_t, Value_t> {
        auto valuePredicate = [](Value_t value) -> Value_t { return value; };
        return ToLookup<decltype(keyPredicate), decltype(valuePredicate), NewKey_t, Value_t>(keyPredicate, valuePredicate);
      }

    };

    template <typename Value_t,
              typename Iterable_t = std::vector<Value_t>>
    auto From(const Value_t* begin, const Value_t* end) -> Enumerable_t<Iterable_t> {
      return Iterable_t(begin, end);
    }

    template <typename Iterator_t,
              typename Iterable_t = std::vector<typename Iterator_t::value_type>>
    auto From(const Iterator_t& begin, const Iterator_t& end) -> Enumerable_t<Iterable_t> {
      return Iterable_t(begin, end);
    }

    template <typename Iterable_t>
    auto From(const Iterable_t& iterable) -> Enumerable_t<Iterable_t> {
      return Iterable_t(iterable.begin(), iterable.end());
    }

  }

}

#endif
