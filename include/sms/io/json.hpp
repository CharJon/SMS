#ifndef SMS_JSON_HPP
#define SMS_JSON_HPP

#include "nlohmann/json.hpp"

namespace sms {

void jsonToFile(const nlohmann::json &j, const std::string &path, int width = 4);

} // namespace sms

#endif // SMS_JSON_HPP
