/*################################################################################
  ##
  ##   Copyright (C) 2016-2018 Keith O'Hara
  ##
  ##   This file is part of the OptimLib C++ library.
  ##
  ##   Licensed under the Apache License, Version 2.0 (the "License");
  ##   you may not use this file except in compliance with the License.
  ##   You may obtain a copy of the License at
  ##
  ##       http://www.apache.org/licenses/LICENSE-2.0
  ##
  ##   Unless required by applicable law or agreed to in writing, software
  ##   distributed under the License is distributed on an "AS IS" BASIS,
  ##   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  ##   See the License for the specific language governing permissions and
  ##   limitations under the License.
  ##
  ################################################################################*/

#ifndef OPTIMLIB_INCLUDES
#define OPTIMLIB_INCLUDES

#include "misc/optim_options.h"

namespace optim
{
    // structs
    #include "misc/optim_structs.h"

    // misc files
    #include "misc/misc.h"

    // line search
    #include "line_search/more_thuente.h"

    // unconstrained optimization
    #include "unconstrained/bfgs.h"
    #include "unconstrained/lbfgs.h"
    #include "unconstrained/newton.h"
    #include "unconstrained/cg.h"
    #include "unconstrained/gd.h"
    #include "unconstrained/de.h"
    #include "unconstrained/de_prmm.h"
    #include "unconstrained/nm.h"
    #include "unconstrained/pso.h"
    #include "unconstrained/pso_dv.h"

    // constrained optimization
    #include "constrained/sumt.h"

    // solving systems of nonlinear equations
    #include "zeros/broyden.h"
}

#endif
