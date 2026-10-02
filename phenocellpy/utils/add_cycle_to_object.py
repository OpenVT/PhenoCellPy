"""
BSD 3-Clause License

Copyright (c) 2023, Juliano Ferrari Gianlupi
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

from .. import phenotypes
from ..types import PhenotypeConfig
from typing import Union


# from .. import phenotypes


def add_cycle_to_object(o: object, phenotype: Union[str, phenotypes.Phenotype], config: PhenotypeConfig = None):
    """
    Attaches a phenotype to `o` as `o.phenotype`.

    :param o: Object to attach the phenotype to. Must support custom attributes
    :param phenotype: An initialized phenotype, or the name of a pre-defined phenotype (see `phenotypes.cycle_names`)
    :type phenotype: str or phenotypes.Phenotype
    :param config: Parameters for the pre-defined phenotype. If None, the phenotype's default config is used. Only
        used if `phenotype` is a name
    :type config: PhenotypeConfig
    """

    if not hasattr(o, "__dict__"):
        raise AttributeError("phenotype class can only be attached to objects that support custom attributes. Object "
                             f"o={o} can't have custom attributes.")

    if not isinstance(phenotype, phenotypes.Phenotype) and phenotype not in phenotypes.cycle_names:
        raise ValueError("Expected `phenotype` parameter to be either an initialized instance of phenotype class or the name "
                         "of a "
                         f"pre-defined phenotype. Got {phenotype}")

    if type(phenotype) == str:
        phenotype = phenotypes.get_phenotype_by_name(phenotype)
        phenotype = phenotype(config)

    setattr(o, "phenotype", phenotype)


def add_phenotype_to_CC3D_cell(cell, phenotype: Union[str, phenotypes.Phenotype], config: PhenotypeConfig = None):
    """
    Attaches a copy of a phenotype to a CompuCell3D cell as `cell.dict["phenotype"]`.

    :param cell: CompuCell3D cell
    :param phenotype: An initialized phenotype, or the name of a pre-defined phenotype (see `phenotypes.cycle_names`)
    :type phenotype: str or phenotypes.Phenotype
    :param config: Parameters for the pre-defined phenotype. If None, the phenotype's default config is used. Only
        used if `phenotype` is a name
    :type config: PhenotypeConfig
    """

    if not hasattr(cell, "dict"):
        raise AttributeError("phenotype class is currently attached to the cell dictionary (i.e., cell.dict), however"
                             f"object {cell} has no dict.")

    if not isinstance(phenotype, phenotypes.Phenotype) and phenotype not in phenotypes.cycle_names:
        raise ValueError("Expected `phenotype` parameter to be either an initialized instance of phenotype class or the name "
                         "of a "
                         f"pre-defined phenotype. Got {phenotype}")

    if type(phenotype) == str:
        phenotype = phenotypes.get_phenotype_by_name(phenotype)
        phenotype = phenotype(config)

    cell.dict["phenotype"] = phenotype.copy()


