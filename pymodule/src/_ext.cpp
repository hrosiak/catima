#include <stdexcept>
#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include "catima/catima.h"
#include "catima/srim.h"
#include "catima/nucdata.h"
#include "catima/convert.h"
#include <iostream>
#include <string>

namespace py = pybind11;
using namespace catima;

std::string catima_info(){
    return "CATIMA version = 1.7\n";
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
    for (size_t i = 0; i < energy_table.size(); i++)
    {
        r.append(energy_table[i]);
    }

    //for(auto e : energy_table){
        //r.append(e);
    //}
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
                    d["sigma_x"] = r.sigma_x;
                    d["cov"] = r.cov;
                    d["tof"] = r.tof;
                    d["sp"] = r.sp;
                    return d;
                    }

PYBIND11_MODULE(_ext,m){
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

     py::class_<Phasespace>(m, "Phasespace")
             .def(py::init<>(),"constructor")
             .def_readwrite("sigma_x", &Phasespace::sigma_x)
             .def_readwrite("sigma_a", &Phasespace::sigma_a)
             .def_readwrite("cov_x", &Phasespace::cov_x);


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
             .def("thickness_cm",py::overload_cast<>(&Material::thickness_cm, py::const_),"get thickness in cm unit")
             .def("thickness_cm",py::overload_cast<double>(&Material::thickness_cm),"set thickness in cm unit")
             .def("I",py::overload_cast<>(&Material::I, py::const_), "get I")
             .def("I",py::overload_cast<double>(&Material::I), "set I")
             .def("number_density",py::overload_cast<>(&Material::number_density, py::const_),"get number density of atoms in cm3 in 10^23 unit")
             .def("number_density",py::overload_cast<int>(&Material::number_density, py::const_),"get number density of atoms of i-th element in cm3 in 10^23 unit")
             .def("number_density_cm2",py::overload_cast<>(&Material::number_density_cm2, py::const_),"get number density of atoms in cm2 in 10^23 unit")
             .def("number_density_cm2",py::overload_cast<int>(&Material::number_density_cm2, py::const_),"get number density of atoms of i-th element in cm2 in 10^23 unit")
	     .def("__str__",&material_to_string);

     py::class_<Layers>(m,"Layers")
             .def(py::init<>(),"constructor")
             .def("add",py::overload_cast<Material>(&Layers::add))
             .def("add_layers",py::overload_cast<const Layers&>(&Layers::add))
             .def("num",&Layers::num)
             .def("thickness",&Layers::thickness)
             .def("thickness_cm",&Layers::thickness_cm)
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
             .def_readwrite("sigma_x", &Result::sigma_x)
             .def_readwrite("cov", &Result::cov)
             .def_readwrite("tof", &Result::tof)
             .def_readwrite("sp", &Result::sp)
             .def("get_dict",&get_result_dict)
             .def("__repr__",[](const Result &self){
                 return py::str(get_result_dict(self));
             });

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
                })
             .def("__repr__",[](const MultiResult &r){
                 py::dict d;
                py::list p;
                d["total_result"] = get_result_dict(r.total_result);
                for(auto& entry:r.results){
                    p.append(get_result_dict(entry));
                    }
                d["results"] = p;
                return py::str(d);
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

    py::enum_<scattering_types>(m,"scattering_types")
            .value("fermi_rossi", scattering_types::fermi_rossi)
            .value("dhighland", scattering_types::dhighland)
            .value("gottschalk", scattering_types::gottschalk)
            .value("atima_scattering", scattering_types::atima_scattering);


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
		 .value("CO2", material::CO2)
		 .value("Methane", material::Methane)
		 .value("Methanol", material::Methanol)
		 .value("Acetone", material::Acetone)
		 .value("Acetylene", material::Acetylene)
		 .value("Adenine", material::Adenine)
		 .value("Adipose_Tissue", material::Adipose_Tissue)
		 .value("Alanine", material::Alanine)
		 .value("Bakelite", material::Bakelite)
		 .value("AgBr", material::AgBr)
		 .value("AgCl", material::AgCl)
		 .value("AgI", material::AgI)
		 .value("Al2O3", material::Al2O3)
		 .value("Amber", material::Amber)
		 .value("Ammonia", material::Ammonia)
		 .value("Aniline", material::Aniline)
		 .value("Anthracene", material::Anthracene)
		 .value("A_150", material::A_150)
		 .value("B_100", material::B_100)
		 .value("BaF2", material::BaF2)
		 .value("BaSO4", material::BaSO4)
		 .value("Benzene", material::Benzene)
		 .value("BeO", material::BeO)
		 .value("BGO", material::BGO)
		 .value("Blood_ICRP", material::Blood_ICRP)
		 .value("Bone_Compact", material::Bone_Compact)
		 .value("Bone_Cortical", material::Bone_Cortical)
		 .value("Brain_ICRP", material::Brain_ICRP)
		 .value("B4C", material::B4C)
		 .value("BC_400", material::BC_400)
		 .value("nButanol", material::nButanol)
		 .value("C_552", material::C_552)
		 .value("CdTe", material::CdTe)
		 .value("CdWO4", material::CdWO4)
		 .value("CaCO3", material::CaCO3)
		 .value("CaF2", material::CaF2)
		 .value("CaO", material::CaO)
		 .value("CaWO4", material::CaWO4)
		 .value("CsF", material::CsF)
		 .value("CsI", material::CsI)
		 .value("CCl4", material::CCl4)
		 .value("Tetrachloroethylene", material::Tetrachloroethylene)
		 .value("Cellophane", material::Cellophane)
		 .value("Chlorobenzene", material::Chlorobenzene)
		 .value("Chloroform", material::Chloroform)
		 .value("Cyclohexane", material::Cyclohexane)
		 .value("Concrete", material::Concrete)
		 .value("Diethyl_Ether", material::Diethyl_Ether)
		 .value("Ethane", material::Ethane)
		 .value("Ethanol", material::Ethanol)
		 .value("Ethylene", material::Ethylene)
		 .value("Eye_lens", material::Eye_lens)
		 .value("Fe2O3", material::Fe2O3)
		 .value("FeO", material::FeO)
		 .value("Freon_12", material::Freon_12)
		 .value("Freon_12B2", material::Freon_12B2)
		 .value("Freon_13", material::Freon_13)
		 .value("Freon_13B1", material::Freon_13B1)
		 .value("Freon_13I1", material::Freon_13I1)
		 .value("Gd2O2S", material::Gd2O2S)
		 .value("GaAs", material::GaAs)
		 .value("Gel_Photo_Emulsion", material::Gel_Photo_Emulsion)
		 .value("Glass_Pyrex", material::Glass_Pyrex)
		 .value("Glass_Lead", material::Glass_Lead)
		 .value("Glucose", material::Glucose)
		 .value("Glutamine", material::Glutamine)
		 .value("Glycerol", material::Glycerol)
		 .value("Guanine", material::Guanine)
		 .value("Gypsum", material::Gypsum)
		 .value("nHeptane", material::nHeptane)
		 .value("nHexane", material::nHexane)
		 .value("KI", material::KI)
		 .value("K2O", material::K2O)
		 .value("LaBr3", material::LaBr3)
		 .value("LaOBr", material::LaOBr)
		 .value("La2O2S", material::La2O2S)
		 .value("Lung", material::Lung)
		 .value("MgCO3", material::MgCO3)
		 .value("MgF2", material::MgF2)
		 .value("MgO", material::MgO)
		 .value("MS20_Tissue", material::MS20_Tissue)
		 .value("Muscle_skeletal", material::Muscle_skeletal)
		 .value("Muscle_strained", material::Muscle_strained)
		 .value("Muscle_sucrose", material::Muscle_sucrose)
		 .value("Muscle_no_sucrose", material::Muscle_no_sucrose)
		 .value("Na2CO3", material::Na2CO3)
		 .value("NaI", material::NaI)
		 .value("NaCl", material::NaCl)
		 .value("Na2O", material::Na2O)
		 .value("NaNO3", material::NaNO3)
		 .value("Naphthalene", material::Naphthalene)
		 .value("Nitrobenzene", material::Nitrobenzene)
		 .value("N2O", material::N2O)
		 .value("Octane", material::Octane)
		 .value("Paraffin", material::Paraffin)
		 .value("nPentane", material::nPentane)
		 .value("PhotoEmulsion", material::PhotoEmulsion)
		 .value("PuO2", material::PuO2)
		 .value("Polyacrylonitrile", material::Polyacrylonitrile)
		 .value("Polycarbonate", material::Polycarbonate)
		 .value("PMMA", material::PMMA)
		 .value("POM", material::POM)
		 .value("Polypropylene", material::Polypropylene)
		 .value("Polystyrene", material::Polystyrene)
		 .value("Propane", material::Propane)
		 .value("nPropanol", material::nPropanol)
		 .value("PVC", material::PVC)
		 .value("Pyridine", material::Pyridine)
		 .value("SiO2", material::SiO2)
		 .value("Skin", material::Skin)
		 .value("Sucrose", material::Sucrose)
		 .value("Teflon", material::Teflon)
		 .value("TlCl", material::TlCl)
		 .value("Toluene", material::Toluene)
		 .value("Trichloroethylene", material::Trichloroethylene)
		 .value("WF6", material::WF6)
		 .value("UC2", material::UC2)
		 .value("UC", material::UC)
		 .value("UO2", material::UO2)
		 .value("Urea", material::Urea)
		 .value("Valine", material::Valine)
		 .value("Iodonaphthalene", material::Iodonaphthalene)
		 .value("C21H24O4", material::C21H24O4)
		 .value("CoRe_Alloy", material::CoRe_Alloy)
		 .value("LLZO_electrolyte", material::LLZO_electrolyte)
		 .value("Nylon", material::Nylon);



    py::class_<Config>(m,"Config")
            .def(py::init<>(),"constructor")
            .def_readwrite("z_effective", &Config::z_effective)
            .def_readwrite("corrections", &Config::corrections)
            .def_readwrite("calculation", &Config::calculation)
            .def_readwrite("low_energy", &Config::low_energy)
            .def_readwrite("scattering", &Config::scattering)
            .def("get",[](const Config &r){
               py::dict d;
               d["z_effective"] = r.z_effective;
               d["corrections"] = r.corrections;
               d["calculation"] = r.calculation;
               d["low_energy"] = r.low_energy;
               d["scattering"] = r.scattering;
               return d;
               })
            .def("__str__",[](const Config &r){
                std::string s;
                s = "z_effective = "+std::to_string(r.z_effective);
                s += ", corrections = "+std::to_string(r.corrections);
                s += ", calculation = "+std::to_string(r.calculation);
                s += ", low_energy = "+std::to_string(r.low_energy);
                s += ", scattering = "+std::to_string(r.scattering);
                return s;
            });

    m.def("angular_scattering_power",py::overload_cast<const Projectile&, const Material&, double>(&angular_scattering_power),"angular scattering power in rad^2/g/cm^2",py::arg("projectile"),py::arg("material"),py::arg("Es2")=Es2_FR);
    m.def("radiation_length",py::overload_cast<const Material&>(radiation_length));
    m.def("srim_dedx_e",&srim_dedx_e);
    m.def("sezi_dedx_e",&sezi_dedx_e, "sezi_dedx_e",  py::arg("projectile"), py::arg("material"), py::arg("config")=default_config);
    m.def("calculate",py::overload_cast<Projectile, const Material&, const Config&>(&calculate),"calculate",py::arg("projectile"), py::arg("material"), py::arg("config")=default_config);
    m.def("calculate",py::overload_cast<const Projectile&, const Layers&, const Config&>(&calculate),"calculate",py::arg("projectile"), py::arg("layers"), py::arg("config")=default_config);
    m.def("calculate",py::overload_cast<const Projectile&, const Phasespace&, const Layers&, const Config&>(&calculate),"calculate",py::arg("projectile"), py::arg("phasespace"),py::arg("layers"), py::arg("config")=default_config);
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
    m.def("save_mocadi", &save_mocadi,py::arg("filename"),py::arg("projectile"),py::arg("layers"),py::arg("psx")=Phasespace(), py::arg("psy")=Phasespace());
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
