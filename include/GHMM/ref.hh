// Copyright (c) 2005 Tobias Sargeant
// 
// Permission is hereby granted, free of charge, to any person obtaining
// a copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject
// to the following conditions:
// 
// The above copyright notice and this permission notice shall be
// included in all copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
// IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR
// ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
// CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#ifndef GHMM_REF_HH_INCLUDED
#define GHMM_REF_HH_INCLUDED

#include <iostream>
#include <stdlib.h>

template<typename T>
struct RefCounter {
  typedef T *ptr_type;
  typedef T ref_type;

  static void incref(T * const &t) { t->incref(); }
  static void decref(T * const &t) { t->decref(); }
  static T &deref(T * const &t) { return *t; }
};

class RefObj {
  mutable int __refcount;
  template<class T> friend class RefCounter;

  RefObj(const RefObj &);
  RefObj &operator=(const RefObj &);

  void incref() const {
    __refcount++;
  }
  void decref() const {
    !--__refcount;
    if (!__refcount) {
      delete this;
    }
  }

protected:
  virtual ~RefObj() {
  }

public:
  int refcount() const {
    return __refcount;
  }
  RefObj() : __refcount(0) {
  }
};

#ifdef __OBJC__
#import <Foundation/Foundation.h>

template<>
struct RefCounter<id> {
  typedef id ptr_type;
  typedef id ref_type;

  static void incref(id const &t) { [t retain]; }
  static void decref(id const &t) { [t release]; }
  static id &deref(id const &t) { return const_cast<id>(t); }
};

#endif

template<typename T, typename R = RefCounter<T> >
class Ref {
  typename R::ptr_type __p;
public:
  Ref() : __p(NULL) {
  }
  template<typename U, typename S>
  Ref(const Ref<U, S> &p) : __p(NULL) {
    *this = p;
  }
  Ref(const Ref &p) : __p(NULL) {
    *this = p;
  }
  Ref(typename R::ptr_type p) : __p(p) {
    if (__p) R::incref(__p);
  }
  ~Ref() {
    if (__p) R::decref(__p);
  }
  template<typename U, typename S>
  Ref &operator=(const Ref<U, S> &p) {
    if (p.ptr()) S::incref(p.ptr());
    if (__p) R::decref(__p);
    __p = p.ptr();
    return *this;
  }
  Ref &operator=(const Ref &p) {
    if (p.__p) R::incref(p.__p);
    if (__p) R::decref(__p);
    __p = p.__p;
    return *this;
  }
  Ref &operator=(typename R::ptr_type p) {
    if (p) R::incref(p);
    if (__p) R::decref(__p);
    __p = p;
    return *this;
  }
  typename R::ref_type &operator*() { return R::deref(__p); }
  const typename R::ref_type &operator*() const { return R::deref(__p); }

  typename R::ptr_type operator->() { return __p; }
  const typename R::ptr_type operator->() const { return __p; }

  typename R::ptr_type ptr() { return __p; }
  const typename R::ptr_type ptr() const { return __p; }

  int refcount() const { return __p->refcount(); }

  template<typename U, typename S>
  bool operator==(const Ref<U, S> &r) const { return r.ptr() == __p; }
  template<typename U, typename S>
  bool operator!=(const Ref<U, S> &r) const { return r.ptr() != __p; }

  bool operator==(const typename R::ptr_type r) const { return r == __p; }
  bool operator!=(const typename R::ptr_type r) const { return r != __p; }
};

#endif
