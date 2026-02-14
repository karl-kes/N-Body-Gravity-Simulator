#pragma once

#include "Body.hpp"
#include "Particle/Particle.hpp"

#include <fstream>
#include <cstdint>
#include <cstring>
#include <vector>
#include <string>
#include <stdexcept>

/*
    Binary format:
    ── Header ──
    uint64_t  num_bodies
    For each body:
        char[32] name  (null-padded)

    ── Per frame ──
    uint64_t  step
    double    time_s
    For each body:
        double x, y, z, vx, vy, vz   (meters, m/s)
*/

class Binary_Output {
private:
    std::ofstream file_;
    std::size_t num_bodies_;
    std::vector<double> buffer_;

public:
    explicit Binary_Output( std::string const &path, Body const *body_list, std::size_t const num )
    : file_{ path, std::ios::binary | std::ios::out }
    , num_bodies_{ num } {
        if ( !file_.is_open() ) {
            throw std::runtime_error{ "Failed to open file: " + path };
        }

        // Header: number of bodies
        uint64_t const n{ static_cast<uint64_t>( num ) };
        file_.write( reinterpret_cast<char const*>( &n ), sizeof( n ) );

        // Header: body names as fixed 32-byte fields
        for ( std::size_t i = 0; i < num; ++i ) {
            char name_buf[32] = {};
            std::strncpy( name_buf, body_list[i].name, 31 );
            file_.write( name_buf, 32 );
        }

        buffer_.reserve( 2 + num * 6 );
    }

    void write(
        Particles const &particles,
        std::size_t const num,
        std::size_t const step,
        double const time ) {

        buffer_.clear();

        // Frame header — store step as reinterpreted double for uniform buffer
        uint64_t const s{ step };
        double step_bits;
        std::memcpy( &step_bits, &s, sizeof( double ) );
        buffer_.push_back( step_bits );
        buffer_.push_back( time );

        // Per-body state
        for ( std::size_t i = 0; i < num; ++i ) {
            buffer_.push_back( particles.pos_x(i) );
            buffer_.push_back( particles.pos_y(i) );
            buffer_.push_back( particles.pos_z(i) );
            buffer_.push_back( particles.vel_x(i) );
            buffer_.push_back( particles.vel_y(i) );
            buffer_.push_back( particles.vel_z(i) );
        }

        file_.write( reinterpret_cast<char const*>( buffer_.data() ),
                     buffer_.size() * sizeof( double ) );
    }

    ~Binary_Output() {
        file_.flush();
    }
};