#pragma once

#include <chrono>
#include <utility>

namespace seal
{
    namespace util
    {
        namespace global_variables
        {
            // inline (not just header-scope): this header is included from many
            // translation units via seal.h -- without inline, each TU gets its own
            // definition of these, and the linker fails with multiple-definition
            // errors as soon as more than one gets linked together.
            inline std::pair<std::chrono::_V2::system_clock::time_point, std::chrono::_V2::system_clock::time_point>
                time_points;

            inline size_t counter;
        } // namespace global_variables

        inline void start_timer()
        {
            global_variables::time_points.first = std::chrono::high_resolution_clock::now();
        }

        inline void stop_timer()
        {
            global_variables::time_points.second = std::chrono::high_resolution_clock::now();
        }

        inline std::chrono::milliseconds get_elapsed_milliseconds()
        {
            return std::chrono::duration_cast<std::chrono::milliseconds>(
                global_variables::time_points.second - global_variables::time_points.first);
        }

        inline void reset_counter()
        {
            global_variables::counter = 0;
        }

        inline void counting()
        {
            global_variables::counter++;
        }

        inline size_t get_count()
        {
            return global_variables::counter;
        }

    } // namespace util
} // namespace seal
