#pragma once

#include <string>

struct StringWrapper
{
  const char * pointer;
  const size_t length;
  StringWrapper() = delete;
  StringWrapper(const StringWrapper&) = delete;
  void operator=(const StringWrapper&) = delete;
  StringWrapper(const std::string& s):
    pointer(s.c_str()),
    length(s.length())
  {}
};
