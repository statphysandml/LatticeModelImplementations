//
// Created by lukas on 05.08.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_LATTICE_UPDATE_HPP
#define LATTICEMODELIMPLEMENTATIONS_LATTICE_UPDATE_HPP


#include "param_helper/params.hpp"


// ToDo: auto update_lattice_site(UpdateFormalism& f, Args&... args) with the &, the neighbours are destroyed after an update -> still true?

// https://www.grimm-jaud.de/index.php/blog/perfect-forwarding
// https://stackoverflow.com/questions/3582001/what-are-the-main-purposes-of-using-stdforward-and-which-problems-it-solves
template<typename UpdateFormalism, typename... Args>
auto update_lattice_site(UpdateFormalism& f, Args&&... args) {
    return f(std::forward<Args>(args)...);
}

template<typename UpdateFormalism, typename... Args>
void global_lattice_update(UpdateFormalism& f, Args&&... args) {
    f(std::forward<Args>(args)...);
}

// ToDo: Introduce something similar for site_update?


class LatticeUpdateParameters : public Parameters {
public:
    using Parameters::Parameters;

    const void write_to_file(const std::string& root_dir) {
        Parameters::write_to_file(root_dir, "lattice_update_params");
    }
};


template<typename Derived>
class LatticeUpdate
{
public:
    template<typename Lattice>
    void initialize (const Lattice &lattice)
    {
        return lattice_update().initialize_update(lattice);
    }

    template<typename Lattice>
    void operator() (Lattice &lattice, uint measure_interval=1)
    {
        return lattice_update().update(lattice, measure_interval);
    }
private:
    Derived& lattice_update() {
        return *static_cast<Derived*>(this);
    }

    const Derived& lattice_update() const {
        return *static_cast<const Derived*>(this);
    }
};

#endif //LATTICEMODELIMPLEMENTATIONS_LATTICE_UPDATE_HPP
