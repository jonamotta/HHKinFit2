#ifndef HHLimitSettingException_
#define HHLimitSettingException_

#include "HHKinFit2/HHKinFit2Core/interface/exceptions/HHEnergyRangeException.h"

#include <iostream>
#include <stdexcept>
#include <string>

namespace HHKinFit2{
  class HHLimitSettingException: public HHEnergyRangeException
{
public:
 HHLimitSettingException(const std::string message)
   : HHEnergyRangeException(message) {};
};
}
#endif
