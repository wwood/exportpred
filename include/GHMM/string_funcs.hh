// Copyright (c) 2005 The Walter and Eliza Hall Institute
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
#ifndef STRING_FUNCS_HH_INCLUDED
#define STRING_FUNCS_HH_INCLUDED

#include <assert.h>
#include <stdarg.h>

#include <string>
#include <vector>

namespace string {
  static inline std::string strip(const std::string &s) {
    int x, y;
    for (x = 0; x < (int)s.size() && isspace(s[x]); x++);
    for (y = (int)s.size() - 1; y >=0 && isspace(s[y]); y--);
    return std::string(s, x, y - x + 1);
  }

  static inline std::vector<std::string> split(const std::string &s, int split = -1, int nsplits = -1) {
    std::vector<std::string> result;
    int x = 0, rx = 0;
    if (split == -1) {
      int y = (int)s.size() - 1;
      while (isspace(s[x])) x++;
      rx = x;
      while (y >= 0 && isspace(s[y])) y--;
      y++;

      while (x < y) {
        if (isspace(s[x])) {
          result.push_back(std::string(s, rx, x - rx));
          while (x < y && isspace(s[x])) x++;
          rx = x;
          if (!--nsplits) {
            result.push_back(std::string(s, rx));
            return result;
          }
        } else {
          x++;
        }
      }
      if (rx != y) {
        result.push_back(std::string(s, rx, y - rx));
      }
    } else {
      while (x < (int)s.size()) {
        if (s[x] == (char)split) {
          result.push_back(std::string(s, rx, x - rx));
          rx = x + 1;
          if (!--nsplits) {
            result.push_back(std::string(s, rx));
            return result;
          }
        }
        x++;
      }
      if (rx < (int)s.size()) {
        result.push_back(std::string(s, rx));
      }
    }
    return result;
  }

  static inline std::string vformat(const char *fmt, va_list args) {
    char *buf;
    if (vasprintf(&buf, fmt, args) == -1) {
      return "";
    } else {
      std::string result = buf;
      free(buf);
      return result;
    }
  }

  static inline std::string format(const char *fmt, ...) {
    va_list args;
    va_start(args, fmt);
    std::string result = vformat(fmt, args);
    va_end(args);
    return result;
  }
};

#endif
