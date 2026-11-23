// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Central definition of the Monte Carlo scalar type.
//
// The CT-INT CPU walker currently supports double precision only.  Keeping the alias here makes
// the precision choice explicit and prevents accidental mixing of float/double in the CT-INT
// implementation.

#ifndef DCA_CONFIG_MC_OPTIONS_HPP
#define DCA_CONFIG_MC_OPTIONS_HPP

namespace dca {
namespace config {

struct McOptions {
  using MCScalar = double;
};

}  // namespace config
}  // namespace dca

#endif  // DCA_CONFIG_MC_OPTIONS_HPP
