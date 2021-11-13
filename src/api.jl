using NLPModels
using LinearAlgebra

abstract type AbstractGNEP{T, S} end

mutable struct GNEP{T, S, NLP <: AbstractNLPModel{T, S}} <: AbstractGNEP{T, S}

 N :: Int64 #number of player
 n :: Array{Int64,1} #array of integer (size of the player strategy)
 player :: Array{NLP, 1} #array of nlpmodels

 x0     :: S
 sol    :: Function #fonction returning a boolean if x is a solution
 class  :: String
 gnsc   :: Bool #true if it is a GNSC
 proj   :: Bool #true if the projection function is furnished
 fnproj :: Function #projection over the feasible set

 function GNEP(N      :: Int64,
               n      :: Array{Int64,1},
               player :: Array{NLP, 1};
               x0     :: S = ones(sum(n)),
               sol    :: Function = x->true,
               class  :: String = "GNEP",
               gnsc   :: Bool   = false,
               proj   :: Bool   = false,
               fnproj :: Function = x -> x) where {NLP, S}

  if length(n)      != N throw(error("wrong size")) end
  if length(x0)     != sum(n) throw(error("wrong size")) end
  if length(player) != N throw(error("wrong size")) end

  return new{eltype(x0), S, NLP}(N, n, player, x0, sol, class, gnsc, proj, fnproj)
 end
end

function pseudo_grad(gnep :: AbstractGNEP{T, S}, x :: S) where {T, S}
 g = zeros(T, sum(gnep.n))
 for i=1:gnep.N
  nI=1+sum(gnep.n[1:i-1]):sum(gnep.n[1:i])
  g[nI] = grad(gnep.player[i],x)[nI]
 end
 return g
end

function pseudo_cons(gnep :: AbstractGNEP{T, S}, x :: S) where {T, S}
 cx = zeros(T, nb_cons(gnep))
 c=0
 for i=1:gnep.N
 ncon = gnep.player[i].meta.ncon
  nI=c+1:c+ncon
  cx[nI] = cons(gnep.player[i],x)
  c+=ncon
 end
 return cx
end

function pseudo_res(gnep :: AbstractGNEP{T, S}, x :: S) where {T, S}
 cx = zeros(T, nb_cons(gnep))
 c=0
 for i=1:gnep.N
 ncon = gnep.player[i].meta.ncon
  nI=c+1:c+ncon
  cx[nI] = cons(gnep.player[i],x)-gnep.player[i].meta.ucon
  c+=ncon
 end
 return cx
end

function pseudo_vio(gnep :: AbstractGNEP{T, S}, x :: S) where {T, S}
 cx = zeros(T, nb_cons(gnep))
 c=0
 for i=1:gnep.N
 ncon = gnep.player[i].meta.ncon
  nI=c+1:c+ncon
  cx[nI] = max.(cons(gnep.player[i],x)-gnep.player[i].meta.ucon,0) + max.(-cons(gnep.player[i],x)+gnep.player[i].meta.lcon,0)
  c+=ncon
 end
 return cx
end

function pseudo_hess(gnep :: AbstractGNEP{T, S}, x :: S) where {T, S}
 H = zeros(T, sum(gnep.n),sum(gnep.n))
 for i=1:gnep.N
  nI=1+sum(gnep.n[1:i-1]):sum(gnep.n[1:i])
  H[nI,nI] = hess(gnep.player[i],x)[nI,nI]
 end
 return H
end

function nb_cons(gnep :: AbstractGNEP)
  count = 0
  for i=1:gnep.N
  count += gnep.player[i].meta.ncon
 end
 return count
end

function pseudo_jac(gnep :: AbstractGNEP{T, S}, x :: S) where {T, S}
 J = zeros(T, nb_cons(gnep),sum(gnep.n))
 c = 0
 for i=1:gnep.N
  nI=1+sum(gnep.n[1:i-1]):sum(gnep.n[1:i])
  ncon = gnep.player[i].meta.ncon
  nC=c+1:c+ncon
  J[nC,nI] = jac(gnep.player[i],x)[:,nI]
  c += ncon
 end
 return J
end

# Lagrange Multiplier of the KKT-GNEP
# (without the bounds constraints)
# (without the signs constraints)
#
function lagrange_multiplier(gnep :: AbstractGNEP{T, S}, x :: S; prec :: Float64 = 1e-6) where {T, S}
 lambda = zeros(T, nb_cons(gnep)+sum(gnep.n))
 c = 0

 for i=1:gnep.N
  nI=1+sum(gnep.n[1:i-1]):sum(gnep.n[1:i])
  nn = length(nI)
  ncon = gnep.player[i].meta.ncon

  cx = cons(gnep.player[i],x)
  ncP = findall(min.(abs.(cx-gnep.player[i].meta.ucon),abs.(gnep.player[i].meta.lcon-cx)) .<= prec)
  ncB = findall(min.(abs.(x[nI]-gnep.player[i].meta.uvar[nI]),abs.(gnep.player[i].meta.lvar[nI]-x[nI])) .<= prec)
  J = vcat(jac(gnep.player[i],x)[ncP,nI],Matrix(1.0I, nn, nn)[ncB,:])
  lambda[c.+union(ncP,ncon.+ncB)] = -LinearAlgebra.pinv(convert(Array,J'))*grad(gnep.player[i],x)[nI]
  c += ncon + nn
 end

 return lambda
end

#bounds problem
#function lagrange_multiplier_bool(gnep :: GNEP, x :: Vector; prec :: Float64 = 1e-6)
# lambda = lagrange_multiplier(gnep, x, prec = prec)
# return (norm(pseudo_jac(gnep,x)'*lambda+pseudo_grad(gnep,x)) <= prec)
#end

function lagrange_multiplier_bool(gnep :: AbstractGNEP{T, S}, x :: S; prec :: Float64 = 1e-6) where {T, S}
 lambda = lagrange_multiplier(gnep, x, prec = prec)
 cl = 0
 res = Inf
 for i=1:gnep.N
  nI=1+sum(gnep.n[1:i-1]):sum(gnep.n[1:i])
  nn = length(nI)
  ncon = gnep.player[i].meta.ncon
  nL=cl+1:cl+ncon+nn
  J = vcat(jac(gnep.player[i],x)[:,nI],Matrix(1.0I, nn, nn))
  res = min(res, norm(J'*lambda[nL]+grad(gnep.player[i],x)[nI]))
  cl += ncon + nn
 end
 return res <= prec
end
