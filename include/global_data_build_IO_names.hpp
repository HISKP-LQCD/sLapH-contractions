#pragma once

#include "global_data.hpp"

/**
 * Fills the random vector and perambulator structs with the paths and file
 * names to read the data.
 */
void build_IO_names(GlobalData &gd, ssize_t const config);
