//
// Created by lukas on 11.10.19.
//

#ifndef MAIN_SITE_MODEL_HPP
#define MAIN_SITE_MODEL_HPP


#include "param_helper/params.hpp"


class SiteModelParameters : public Parameters {
public:
    using Parameters::Parameters;

    void write_to_file(const std::string& root_dir) {
        Parameters::write_to_file(root_dir, "model_params");
    }
};


template <typename ModelCL>
class SiteModel
{
private:
    ModelCL& site_model() {
        return *static_cast<ModelCL*>(this);
    }
};

#endif //MAIN_SITE_MODEL_HPP
