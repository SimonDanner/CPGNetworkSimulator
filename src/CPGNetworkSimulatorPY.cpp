/*  CPGNetworkSimulatorPY.cpp: pybind11 bindings  
    Copyright (C) 2019  Simon Danner

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "CPGNetworkSimulator.hpp"

namespace py = pybind11;

PYBIND11_PLUGIN(_CPGNetworkSimulator) {
    py::module m("_CPGNetworkSimulator", "Network Simulator");
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
        .def("controlled_step", &CPGNetworkSimulator::controlled_step)
        .def("getAct", &CPGNetworkSimulator::getAct)
        .def("setAlpha", &CPGNetworkSimulator::setAlpha)
        .def("updateVariable", &CPGNetworkSimulator::updateVariable)
        .def("setBodyTilt",&CPGNetworkSimulator::setBodyTilt)
        .def("setupVariableVector",&CPGNetworkSimulator::setupVariableVector)
        .def("updateVariableVector",&CPGNetworkSimulator::updateVariableVector)
        .def("updateParameter",&CPGNetworkSimulator::updateParameter)
        .def("setLscond", &CPGNetworkSimulator::setLscond)
        .def("getState", &CPGNetworkSimulator::getState)
        .def("setState", &CPGNetworkSimulator::setState);
    return m.ptr();
}

