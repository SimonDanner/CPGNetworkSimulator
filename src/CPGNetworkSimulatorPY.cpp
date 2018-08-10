#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "CPGNetworkSimulator.hpp"

namespace py = pybind11;

PYBIND11_MODULE(CPGNetworkSimulator, m) {
    py::class_<LimbSensorCondition>(m, "LimbSensorCondition")
        .def(py::init<>())
        .def(py::init<int>())
        .def_readwrite("Ia", &LimbSensorCondition::Ia)
        .def_readwrite("Ib", &LimbSensorCondition::Ib)
        .def_readwrite("II", &LimbSensorCondition::II)
        .def_readwrite("cutaneous", &LimbSensorCondition::cutaneous);

    py::class_<CPGNetworkSimulator>(m, "CPGNetworkSimulator")
        .def(py::init<const std::string, const std::vector<std::string>, const std::vector<std::vector<std::string>>>())
        .def("step", &CPGNetworkSimulator::step)
        .def("getAct", &CPGNetworkSimulator::getAct)
        .def("setAlpha", &CPGNetworkSimulator::setAlpha)
        .def("updateVariable", &CPGNetworkSimulator::updateVariable)
        .def("setBodyTilt",&CPGNetworkSimulator::setBodyTilt)
        .def("setLscond", &CPGNetworkSimulator::setLscond);
}

