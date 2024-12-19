#include "sms/io/json.hpp"

#include <fstream>

namespace sms {

void jsonToFile(const nlohmann::json &j, const std::string &path, int width) {
    std::ofstream file(path);
    file << std::setw(width) << j << std::endl;
}

} // namespace sms