#include <stdexcept>
#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include "catima/catima.h"
#include "catima/srim.h"
#include "catima/nucdata.h"
#include <iostream>
#include <string>
namespace py = pybind11;
using namespace catima;

void catima_info(){
    printf("CATIMA version = 1.5\n");
    printf("number of energy points = %d\n",max_datapoints);
    printf("min energy point = 10^%lf MeV/u\n",logEmin);
    printf("max energy point = 10^%lf MeV/u\n",logEmax);
}

std::string  material_to_string(const Material &r){
                std::string s;
                auto n = r.ncomponents();
                for(int i = 0; i < n; i++){
                    auto el = r.get_element(i);
                    s += "#"+std::to_string(i);
                    s += ": A = "+std::to_string(el.A) + ", Z = "+std::to_string(el.Z)+ ", stn = "+std::to_string(el.stn)+"\n";
                 }
                return s;
            };

py::list storage_info(){
    py::list res;
    for(int i=0; i<max_storage_data;i++){
        auto& data = _storage.Get(i);
        if(data.p.A>0 && data.p.Z && data.m.ncomponents()>0){
            py::list mat;
            py::dict d;
            py::list p;
            p.append(data.p.A);
            p.append(data.p.Z);
            d["projectile"] = p;
            d["matter"] = material_to_string(data.m);
            d["config"] = py::cast(data.config);
            res.append(d);
        }
    }
    return res;
}

py::list get_energy_table(){
    py::list r;
    for(auto e : energy_table){
        r.append(e);
    }
    return r;
}

py::list get_data(Projectile& p, const Material &m, const Config& c=default_config){
    py::list r;
    auto& data = _storage.Get(p, m, c);
    py::list ran;
    py::list rans;
    py::list av;
    for(double e:data.range){ran.append(e);}
    for(double e:data.range_straggling)rans.append(e);
    for(double e:data.angular_variance)av.append(e);
    r.append(ran);
    r.append(rans);
    r.append(av);
    return r;
}

Material py_make_material(py::list d, double density=0.0, double thickness=0.0, double ipot=0.0, double mass=0.0){
    Material m;
    if(density>0.0)m.density(density);
    if(ipot>0.0)m.I(ipot);
    if(mass>0.0)m.M(mass);
    if(thickness>0.0)m.thickness(thickness);
    for(int i=0;i<d.size();i++){
        py::list e(d[i]);
        if(e.size() != 3)throw std::invalid_argument("invalid Material constructor argument");
        double a = e[0].cast<double>();
        int z = e[1].cast<int>();
        double stn = e[2].cast<double>();
        m.add_element(a,z,stn);
    }

    return m;
}

py::dict get_result_dict(const Result& r){
                    py::dict d;
                    d["Ein"] = r.Ein;
                    d["Eout"] = r.Eout;
                    d["Eloss"] = r.Eloss;
                    d["range"] = r.range;
                    d["dEdxi"] = r.dEdxi;
                    d["dEdxo"] = r.dEdxo;
                    d["sigma_E"] = r.sigma_E;
                    d["sigma_r"] = r.sigma_r;
                    d["sigma_a"] = r.sigma_a;
                    d["tof"] = r.tof;
                    d["sp"] = r.sp;
                    return d;
                    }

PYBIND11_MODULE(pycatima,m){
     py::class_<Projectile>(m,"Projectile")
             .def(py::init<>(),"constructor")
             .def(py::init<double, double, double, double>(), "constructor", py::arg("A"),py::arg("Z"),py::arg("Q")=0, py::arg("T")=0)
             .def("__call__",&Projectile::operator())
             .def("A",[](const Projectile& p){return p.A;})
             .def("Z",[](const Projectile& p){return p.Z;})
             .def("Q",[](const Projectile& p){return p.Q;})
             .def("T",[](const Projectile& p){return p.T;})
             .def("T",[](Projectile& p, double v){p.T = v;});
             //.def_readwrite("A", &Projectile::A)
             //.def_readwrite("Z", &Projectile::Z)
             //.def_readwrite("T", &Projectile::T)
             //.def_readwrite("Q", &Projectile::Q);

     py::class_<Target>(m,"Target")
             .def(py::init<>(),"constructor")
             .def_readwrite("A",&Target::A)
             .def_readwrite("Z",&Target::Z)
             .def_readwrite("stn",&Target::stn);


     py::class_<Material>(m,"Material")
             .def(py::init<>(),"constructor")
             .def(py::init<const Material&>(),"constructor")
             .def(py::init<double, int, double, double, double>(),"constructor", py::arg("A"),py::arg("Z"),py::arg("density")=0.0,py::arg("thickness")=0.0,py::arg("i_potential")=0.0)
             .def(py::init(&py_make_material),"constructor", py::arg("elements"),py::arg("density")=0.0,py::arg("thickness")=0.0,py::arg("i_potential")=0.0, py::arg("mass")=0.0)
             .def("add_element",&Material::add_element)
             .def("ncomponents",&Material::ncomponents)
             .def("density",py::overload_cast<>(&Material::density, py::const_), "get density")
             .def("density",py::overload_cast<double>(&Material::density), "set density")
             .def("molar_mass",py::overload_cast<>(&Material::M, py::const_), "get mass")
             .def("thickness",py::overload_cast<>(&Material::thickness, py::const_), "get thickness")
             .def("thickness",py::overload_cast<double>(&Material::thickness), "set thickness")
             .def("thickness_cm",py::overload_cast<double>(&Material::thickness_cm),"set thickness in cm unit")
             .def("I",py::overload_cast<>(&Material::I, py::const_), "get I")
             .def("I",py::overload_cast<double>(&Material::I), "set I")
             .def("__str__",&material_to_string);

     py::class_<Layers>(m,"Layers")
             .def(py::init<>(),"constructor")
             .def("add",py::overload_cast<Material>(&Layers::add))
             .def("add_layers",py::overload_cast<const Layers&>(&Layers::add))
             .def("num",&Layers::num)
//             .def("__getitem__",&Layers::operator[],  py::is_operator())
             .def("__getitem__",[](Layers &r, int i)->Material*
                {
                 if(i>=r.num()){
                     throw std::invalid_argument("index out of range");}
                 return &r[i];
                },  py::is_operator(),py::return_value_policy::automatic_reference)
             .def("get",&Layers::operator[])
             .def(py::self + py::self)
             .def("__add__",[](const Layers s, const Material& m){return s+m;});

     py::class_<Result>(m,"Result")
             .def(py::init<>(),"constructor")
             .def_readwrite("Ein", &Result::Ein)
             .def_readwrite("Eout", &Result::Eout)
             .def_readwrite("Eloss", &Result::Eloss)
             .def_readwrite("range", &Result::range)
             .def_readwrite("dEdxi", &Result::dEdxi)
             .def_readwrite("dEdxo", &Result::dEdxo)
             .def_readwrite("sigma_E", &Result::sigma_E)
             .def_readwrite("sigma_a", &Result::sigma_a)
             .def_readwrite("sigma_r", &Result::sigma_r)
             .def_readwrite("tof", &Result::tof)
             .def_readwrite("sp", &Result::sp)
             .def("get_dict",&get_result_dict);

     py::class_<MultiResult>(m,"MultiResult")
             .def(py::init<>(),"constructor")
             .def_readwrite("total_result", &MultiResult::total_result)
             .def_readwrite("results", &MultiResult::results)
//             .def_readwrite("Eout",&MultiResult::total_result.Eout)
             .def("__getitem__",[](MultiResult &r, int i){
                return py::cast(r.results[i]);
             },py::is_operator())
             .def("__getattr__",[](MultiResult &r, std::string& k){
                 if(k.compare("Eout")==0){
                     return py::cast(r.total_result.Eout);
                 }
                 else if(k.compare("sigma_a")==0){
                     return py::cast(r.total_result.sigma_a);
                 }
                 else if(k.compare("tof")==0){
                     return py::cast(r.total_result.tof);
                 }
                 else if(k.compare("Eloss")==0){
                     return py::cast(r.total_result.Eloss);
                 }
                 else{
                     return py::cast(NULL);
                 }
             },py::is_operator())
             .def("get_dict",[](const MultiResult &r){
                py::dict d;
                py::list p;
                d["result"] = get_result_dict(r.total_result);
                for(auto& entry:r.results){
                    p.append(get_result_dict(entry));
                    }
                d["partial"] = p;
                return d;
                });

    py::enum_<z_eff_type>(m,"z_eff_type")
            .value("none", z_eff_type::none)
            .value("pierce_blann", z_eff_type::pierce_blann)
            .value("anthony_landorf", z_eff_type::anthony_landorf)
            .value("hubert", z_eff_type::hubert)
            .value("winger", z_eff_type::winger)
            .value("schiwietz", z_eff_type::schiwietz)
            .value("cglobal", z_eff_type::global)
            .value("atima14", z_eff_type::atima14);

    py::enum_<corrections>(m,"corrections")
            .value("no_barkas", corrections::no_barkas)
            .value("no_lindhard", corrections::no_lindhard)
            .value("no_shell_correction", corrections::no_shell_correction)
            .value("no_highenergy", corrections::no_highenergy);

    py::enum_<omega_types>(m,"omega_types")
            .value("atima", omega_types::atima)
            .value("bohr", omega_types::bohr);

    py::enum_<low_energy_types>(m,"low_energy_types")
            .value("srim_85", low_energy_types::srim_85)
            .value("srim_95", low_energy_types::srim_95);

    py::enum_<material>(m,"material")
            .value("Plastics", material::Plastics)
            .value("Air", material::Air)
            .value("CH2", material::CH2)
            .value("lH2", material::lH2)
            .value("lD2", material::lD2)
            .value("Water", material::Water)
            .value("Diamond", material::Diamond)
            .value("Glass", material::Glass)
            .value("ALMG3", material::ALMG3)
            .value("ArCO2_30", material::ArCO2_30)
            .value("CF4", material::CF4)
            .value("Isobutane", material::Isobutane)
            .value("Kapton", material::Kapton)
            .value("Mylar", material::Mylar)
            .value("NaF", material::NaF)
            .value("P10", material::P10)
            .value("Polyolefin", material::Polyolefin)
            .value("CmO2", material::CmO2)
            .value("Suprasil", material::Suprasil)
            .value("HAVAR", material::HAVAR)
            .value("Steel", material::Steel)
            .value("CO2", material::CO2);


    py::class_<Config>(m,"Config")
            .def(py::init<>(),"constructor")
            .def_readwrite("z_effective", &Config::z_effective)
            .def_readwrite("corrections", &Config::corrections)
            .def_readwrite("calculation", &Config::calculation)
            .def("get",[](const Config &r){
               py::dict d;
               d["z_effective"] = r.z_effective;
               d["corrections"] = r.corrections;
               d["calculation"] = r.calculation;
               return d;
               })
            .def("__str__",[](const Config &r){
                std::string s;
                s = "z_effective = "+std::to_string(r.z_effective);
                s += ", corrections = "+std::to_string(r.corrections);
                s += ", calculation = "+std::to_string(r.calculation);
                return s;
            });

    m.def("srim_dedx_e",&srim_dedx_e);
    m.def("sezi_dedx_e",&sezi_dedx_e, "sezi_dedx_e",  py::arg("projectile"), py::arg("material"), py::arg("config")=default_config);
    m.def("calculate",py::overload_cast<Projectile, const Material&, const Config&>(&calculate),"calculate",py::arg("projectile"), py::arg("material"), py::arg("config")=default_config);
    m.def("calculate_layers",py::overload_cast<const Projectile&, const Layers&, const Config&>(&calculate),"calculate_layers",py::arg("projectile"), py::arg("material"), py::arg("config")=default_config);
    m.def("dedx_from_range",py::overload_cast<const Projectile&, const Material&, const Config&>(&dedx_from_range),"calculate",py::arg("projectile") ,py::arg("material"), py::arg("config")=default_config);
    m.def("dedx_from_range",py::overload_cast<const Projectile&, const std::vector<double>&, const Material&, const Config&>(&dedx_from_range),"calculate",py::arg("projectile"), py::arg("energy") ,py::arg("material"), py::arg("config")=default_config);
    m.def("dedx",py::overload_cast<const Projectile&, const Material&, const Config&>(&dedx), "dedx",py::arg("projectile"), py::arg("material"), py::arg("config")=default_config);
    m.def("range",py::overload_cast<const Projectile&, const Material&, const Config&>(&range), "range",py::arg("projectile"), py::arg("material"), py::arg("config")=default_config);
    m.def("energy_out",py::overload_cast<const Projectile&, const std::vector<double>&, const Material&, const Config&>(&energy_out),"energy_out",py::arg("projectile"), py::arg("energy") ,py::arg("material"), py::arg("config")=default_config);
    m.def("energy_out",py::overload_cast<const Projectile&, const Material&, const Config&>(&energy_out),"energy_out",py::arg("projectile"), py::arg("material"), py::arg("config")=default_config);
    m.def("lindhard",&bethek_lindhard);
    m.def("lindhard_X",&bethek_lindhard_X);
    m.def("get_material",py::overload_cast<int>(&get_material));
    m.def("get_data",py::overload_cast<Projectile&, const Material&, const Config&>(get_data),"list of data",py::arg("projectile"),py::arg("material"),py::arg("config")=default_config);
    m.def("w_magnification",[](Projectile& p, double energy, const Material& m, const Config& c){
        py::list l;
        auto r = w_magnification(p, energy, m, c);
        l.append(r.first);
        l.append(r.second);
        return l;
    });
    m.def("catima_info",&catima_info);
    m.def("storage_info",&storage_info);
    m.def("get_energy_table",&get_energy_table);
    m.def("energy_table",[](int i){return energy_table(i);});
    m.def("z_effective",&z_effective);
    m.attr("max_datapoints") = max_datapoints;
    m.attr("max_storage_data") = max_storage_data;
    m.attr("logEmin")=logEmin;
    m.attr("logEmax")=logEmax;
}
