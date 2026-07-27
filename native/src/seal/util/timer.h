#pragma once

#include <chrono>
#include <utility>

namespace seal
{
    namespace util
    {
        namespace global_variables
        {
            std::pair<std::chrono::_V2::system_clock::time_point, std::chrono::_V2::system_clock::time_point>
                time_points;

            size_t counter;
        } // namespace global_variables

        void start_timer()
        {
            global_variables::time_points.first = std::chrono::high_resolution_clock::now();
        }

        void stop_timer()
        {
            global_variables::time_points.second = std::chrono::high_resolution_clock::now();
        }

        std::chrono::milliseconds get_elapsed_milliseconds()
        {
            return std::chrono::duration_cast<std::chrono::milliseconds>(
                global_variables::time_points.second - global_variables::time_points.first);
        }

        void reset_counter()
        {
            global_variables::counter = 0;
        }

        void counting()
        {
            global_variables::counter++;
        }

        size_t get_count()
        {
            return global_variables::counter;
        }

    } // namespace util
} // namespace seal