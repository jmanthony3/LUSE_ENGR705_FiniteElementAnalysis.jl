module LUSE_ENGR705_FiniteElementAnalysis

export FEA
export transformationmatrix
export transformationmatrixstar
export prepareelements
export globalstiffnessmatrix
export reducedglobalstiffnessmatrix
export prepareelements_line
export prepareelements_beam
export prepareelements_triangular
export solve
export axialstresses
export shearstresses
include("functions.jl")

end