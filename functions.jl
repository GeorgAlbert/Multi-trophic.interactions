############################################################################
# Collection of functions to simulate dynamics in trophic networks used in #
# Albert et al. (202x): The hidden role of multi-trophic interactions in   #
# driving diversity-productivity relationships. Ecology Letters            #
############################################################################

## Wrapper function to flexibly generate set of starting parameters and starting densities
# Can be directly used to create a ODE problem
# Required arguments are the number of animals ('nA') and of producers ('nP')
# To simulate resource-use dissimilarity, the number of compartments ('nC') needs to be defined and a resource-use dissimilarity scenario needs to be provided (use RUD-function)
# Is flexible to use predefined body masses by providing optional arguments with vector of body masses for animals ('bodymass_animal'), producers ('bodymass_producer') or both ('bodymass')
function startparam(nA, nP; nC = 1, RUDscenario = undef, bodymass = undef, bodymass_producer = undef, bodymass_animal = undef)

	if nA == 1 & nP > 5
		@warn("Finding a suitable food web topology with nA = 1 can take very long if nP is too high!")
	end

	## number of resources
	nN = 2

	## Bodymasses
	if bodymass == undef

		if bodymass_producer == undef
			bodymass_producer_fixed = false
			bodymass_producer = 10 .^ rand(Distributions.Uniform(log(10, 10^0), log(10, 10^6)), nP)
		else
			if (sum(bodymass_producer .< 10^log(10, 10^0)) > 0) | (sum(bodymass_producer .> 10^log(10, 10^6)) > 0)
				@warn("(Some) body masses of producers are outside of model specific range")
			end
			bodymass_producer_fixed = true
		end

		if bodymass_animal == undef
			bodymass_animal_fixed = false
			bodymass_animal = 10 .^ rand(Distributions.Uniform(log(10, 10^2), log(10, 10^12)), nA)
		else
			if (sum(bodymass_animal .< 10^log(10, 10^2)) > 0) | (sum(bodymass_animal .> 10^log(10, 10^12)) > 0)
				@warn("(Some) body masses of animals are outside of model specific range")
			end
			bodymass_animal_fixed = true
		end

		bodymass = [bodymass_animal..., bodymass_producer...]


		# make sure there are no unconnected species OR animals are basal (i.e. don't have resource)
		# define optimal body mass ratio and gamma to check structure
		optimalmassratio = 100
		gamma = 2

		if nA > 0
			if ((bodymass_producer_fixed == false) | (bodymass_animal_fixed == false))
				while checkstruct(bodymass, nP, optimalmassratio, gamma)

					if bodymass_producer_fixed == false
						bodymass_producer = 10 .^ rand(Distributions.Uniform(log(10, 10^0), log(10, 10^6)), nP)
					end
					if bodymass_animal_fixed == false
						bodymass_animal = 10 .^ rand(Distributions.Uniform(log(10, 10^2), log(10, 10^12)), nA)
					end
					bodymass = [bodymass_animal..., bodymass_producer...]
				end
			end
		end
	end

	if length(bodymass) != (nA + nP)
		@warn("Number of provided body masses is different from nA+nP")
	end
	if nA > 0
		if checkstruct(bodymass, nP, optimalmassratio, gamma)
			@warn("Food web topology with isolated species or basal animals")
		end
	end


	## Starting densities:
	# Animals:
	A = rand(Distributions.Uniform(0, 10), nA)
	while sum(A .== 0) > 0
		A[A .== 0] .= rand(Distributions.Uniform(0, 10), sum(A .== 0))
	end
	# producers:
	P = rand(Distributions.Uniform(0, 10), nP)
	while sum(P .== 0) > 0
		P[P .== 0] .= rand(Distributions.Uniform(0, 10), sum(P .== 0))
	end

	# Resource:
	# - Supply concentration
	S = [50, 25]

	N = repeat(map(x -> rand(Distributions.Uniform(S[x]/2, S[x])), 1:nN) ./ nC, inner = nC) # resources are evenly split between compartments



	## remaining parameters
	# maximum specific growth rate (producer)
	r = collect((bodymass[nA+1:end].^-0.25)')

	# content of resource in producer species
	c = [repeat([1], outer = nP) repeat([0.5], outer = nP)]

	# half-saturation constants
	K = reshape(rand(Distributions.Uniform(0.1, 0.2), nP * nN), nN, :)

	# System turnover
	D = 0.25

	# relative consumption rate
	omega = 1 ./ preyspec(bodymass, nP, optimalmassratio, gamma)
	omega[omega .== Inf] .= 0

	# capture coefficient
	b = capturecoef(bodymass, nP, optimalmassratio, gamma)

	# time lost due to consumer interference
	int = sigmalimit(0.8, 0.2; n = length(bodymass))

	# handling time
	h = handling(bodymass)

	# conversion efficiency - herbivory
	eP = 0.545

	# conversion efficiency - carnivory
	eA =  0.906

	# metabolic demand
	x = vcat(0.141 .* bodymass[1:nA].^-0.305, 0.138 .* bodymass[nA+1:end].^-0.25)

	# RUD / compartments
	# - check if (more or less) correctly initialized
	if (nC != 1)
		if (RUDscenario == undef)
			@warn("Need to specify RUDscenario!")
		elseif (size(RUDscenario) != (nP, nC))
			@warn("RUDscenario has wrong dimensions. Check if nP and nC are correct!")
		end
	end
	# if RUDscenario isn't defined we assume nC is 1 and create an array where all producers access that compartment
	RUDscenario == undef && (RUDscenario = Array{Bool}(fill(1, nP)))

	# Hill exponent
	q = hcat(map(x -> (x ./ bodymass), bodymass)...) |>
		x -> 1 .* x.^2 ./ ((10^2)^2 .+ x.^2)



	## Output:
	# Parametervector
	# - indices:
	#	 1  2  3  4  5  6      7  8  9    10 11        12  13  14 15  16  17  18  19          20                21
	# - output vector
	p = (r, S, c, K, D, omega, b, q, int, h, bodymass, eP, eA, x, nA, nP, nN, nC, RUDscenario, optimalmassratio, gamma)
	# - name vector to help identify parameters
	p_string = "p = (r, S, c, K, D, omega, b, q, int, h, bodymass, eP, eA, x, nA, nP, nN, nC, RUDscenario, optimalmassratio, gamma)"
	# starting conditions
	u0 = collect([A' P' N'])

	return [p, u0, p_string]
end



## modify startparameters
# - 'startparam_orig' is the startparmeters to be modified (e.g. output of startparam-function)
# - 'selectP' is index-vector of the plant species to keep
# - 'nA_new' specifies the animal species richness of a set of animals species that are newly generated
# - 'RUD' can be used to introduce new RUD-scenario
function modify_startparam(startparam_orig; selectP = collect(1:startparam_orig[1][16]), nA_new = undef, RUD = undef)

	startparam_copy = copy(startparam_orig)

	if nA_new != undef
		# create startparameters with new animal richness and the same producer body masses
		startparam_tmp = startparam(nA_new,
			startparam_copy[1][16],
			bodymass_producer = startparam_copy[1][11][startparam_copy[1][15]+1:end])

		# replace parameters in 'startparam' relevant for the interaction with new animal species
		startparam_copy[1] = (startparam_copy[1][1:6]...,
			startparam_tmp[1][7:11]...,
			startparam_copy[1][12:13]...,
			startparam_tmp[1][14:15]...,
			startparam_copy[1][16:end]...)
	end


	nP = length(selectP)

	## remove animal species that go extint as they loose resource species
	# create initial indexvector of all animal and selected producer species
	indexvec = collect(1:length(startparam_copy[1][11])) |>
		x -> deleteat!(x, deleteat!(collect(1:startparam_copy[1][16]), selectP) .+ startparam_copy[1][15])

	# remove animal species from index vector that lose resource species
	while sum(preyspec(startparam_copy[1][11][indexvec], nP, startparam_copy[1][20], startparam_copy[1][21]) .> 0) != length(startparam_copy[1][11][indexvec])-nP

		indexvec = [indexvec[preyspec(startparam_copy[1][11][indexvec], nP, startparam_copy[1][20], startparam_copy[1][21]) .!= 0]..., indexvec[end-nP+1:end]...]
	end

	nA = length(indexvec) - nP


	## create new parameter vector
	startparam_copy[1] = (startparam_copy[1][1][:, selectP],
			   startparam_copy[1][2],
			   startparam_copy[1][3][selectP, :],
			   startparam_copy[1][4][:, selectP],
			   startparam_copy[1][5],
			   startparam_copy[1][6], # omega is updated later
			   startparam_copy[1][7][indexvec, indexvec],
			   startparam_copy[1][8][indexvec, indexvec],
			   startparam_copy[1][9][indexvec],
			   startparam_copy[1][10][indexvec, indexvec],
			   startparam_copy[1][11][indexvec],
			   startparam_copy[1][12],
			   startparam_copy[1][13],
			   startparam_copy[1][14][indexvec],
			   nA,
			   nP,
			   startparam_copy[1][17],
			   startparam_copy[1][18],
			   RUD == undef ? startparam_copy[1][19][selectP, :] : RUD[selectP, :],
			   startparam_copy[1][20],
			   startparam_copy[1][21])

	omega_tmp = 1 ./ preyspec(startparam_copy[1][11], nP, startparam_copy[1][20], startparam_copy[1][21])
	omega_tmp[omega_tmp .== Inf] .= 0

	startparam_copy[1] = (startparam_copy[1][1:5]..., omega_tmp, startparam_copy[1][7:end]...)


	## starting densities
	startparam_copy[2] = [startparam_copy[2][indexvec]..., startparam_copy[2][startparam_orig[1][15]+startparam_orig[1][16]+1:end]...]


	return startparam_copy
end



## Define ODE problem
function foodwebsim(du, u, p, t)
	# set densities below threshold to zero
	u[u .<= 10^-6] .= zeros()

	# sum up resources accessible by producer species
	sumnut = hcat(map(y -> map(x -> sum(u[(p[15] + p[16] + 1):end][
	                                    collect((((x - 1) * p[18]) + 1):(x * p[18]))[p[19][y, :]]]),
								1:p[17]),
					  1:p[16])...)

	# calculate growth rates of producers
	grow = growth(sumnut, p[1], p[4])

	# calculate feeding rates
	feed = feedingrate(u[1:(p[15] + p[16])], p[6], p[7], p[8], p[9], p[10], p[11])

	## ODE:
	# animals
	du[1:p[15]] =
		p[12] .* u[1:p[15]] .* map(x -> sum(feed[(end - p[16] + 1):end, 1:p[15]][:, x]), 1:p[15]) .+
		p[13] .* u[1:p[15]] .* map(x -> sum(feed[1:(end - p[16] + 1), 1:p[15]][:, x]), 1:p[15]) .-
		map(x -> sum(u[1:p[15]] .* feed[1:p[15], 1:p[15]][x, :]), 1:p[15]) .-
		p[14][1:p[15]] .* u[1:p[15]]

	# producers
	du[(p[15] + 1):(p[15] + p[16])] =
		grow .* u[(p[15] + 1):(p[15] + p[16])] .-
		map(x -> sum(u[1:p[15]] .* feed[(p[15] + 1):end, 1:p[15]][x, :]), 1:p[16]) .-
		p[14][(p[15] + 1):(p[15] + p[16])] .* u[(p[15] + 1):(p[15] + p[16])]

	# resources
	du[(p[15] + p[16] + 1):end] =
		p[5] .* (repeat(p[2] ./ p[18], inner = p[18]) .- u[(p[15] + p[16] + 1):end]) .-
		map(z -> sum(hcat(map(y ->
							map(x -> p[3][:, x] .* grow .* u[(p[15] + 1):(p[15] + p[16])],
		 					  1:p[17])[y] .*
				 		  	map(x -> p[19] .* u[(p[15] + p[16] + 1):end][((x-1)*p[18] + 1):(x*p[18])]' ./ sumnut'[:, x],
						  	  1:p[17])[y],
							1:p[17])...)[:, z]),
		  1:p[17]*p[18])

	# set density changes below threshold to zero
	du[u .<= 10^-6] .= zeros()
end



## Create scenarios of RUD gradient or random RUD scenario
function RUD(nP, nC; random = false)
	if random
		RUD = reshape(rand(Bool, nP*nC), nP, nC)

		# make sure each producer has at least access to one compartment
		while (sum(maximum(RUD, dims = 1)) + sum(maximum(RUD, dims = 2))) < nP + nC
			RUD = reshape(rand(Bool, nP*nC), nP, nC)
		end
	else
		if !((nP % nC == 0) | (nC % nP == 0))
			@warn("nC or nP is not multiple of the other. Unbalanced RUD scenario!")
		end

		RUD = map(x -> vcat([[repeat([true], x)..., repeat([false], nC-x)...][[collect(nC-i+1:nC)... collect(1:nC-i)...]] for i in (collect(1:nP).%nC).+1]...), 1:nC)

		# reorders, so that producer 1 starts in compartment 1 (just aesthetical fix)
		RUD .= map(x -> x[:, [collect(3:size(x)[2])..., 1, 2]], RUD)
	end

	return RUD
end



## Functions used within previous
# draw values from normal distribution and redraw values outside of range of 3 sigma
function sigmalimit(mu, sigma; n = 1)
	x = rand(Distributions.Normal(mu, sigma), n)

	while sum((x .< mu - 3 * sigma) .| (x .> mu + 3 * sigma)) > 0
		tmp = ((x .< mu - 3 * sigma) .| (x .> mu + 3 * sigma))

		x[tmp] = rand(Distributions.Normal(mu, sigma), sum(tmp))
	end

	return x
end


# check if all have at least one predator or one prey (no isolation) OR if any animal species is basal
# false --> conditions met!
function checkstruct(bodymass, nP, optimalmassratio, gamma)

	if length(bodymass) != nP
		(sum(preyspec(bodymass, nP, optimalmassratio, gamma) + predspec(bodymass, nP, optimalmassratio, gamma) .== zeros()) > 0)     |     (sum(map(x -> sum(feedingeff(bodymass, nP, optimalmassratio, gamma)[:, x]), 1:Int(length(bodymass)-nP)) .== zeros()) > 0)
	else
		missing
	end
end


# calculate producer growth rates
function growth(N, r, K)
  map(x -> minimum(((r .* N) ./ (K .+ N))[:, x]), 1:size(K, 2))
end


# calculate feeding efficiency
function feedingeff(bodymass, nP, optimalmassratio, gamma)

	# make sure bodymass is vector
	if size(bodymass)[1] == 1
		bodymass = bodymass'
	end

	feedingeff = hcat(map(x -> (x ./ bodymass) ./ optimalmassratio, bodymass)...) |>
		x -> (x.*exp.(1 .-x)).^gamma


	feedingeff[feedingeff .< 0.01] .= zeros() # set values below threshold to zero
	feedingeff[:, Int(length(bodymass)-nP+1):end] .= zeros() # make sure producers don't consume other species!

	return feedingeff
end


# calculate number of prey species
function preyspec(bodymass, nP, optimalmassratio, gamma)
	prey = Vector{Float64}()

	for i in 1:length(bodymass)
		tmp = length(findall(feedingeff(bodymass, nP, optimalmassratio, gamma)[:, i] .!= zeros()[]))
		push!(prey, tmp)
	end

	# make sure producers have no prey (highly unlikely)
	prey[Int(end-nP+1):end] .= zeros()

	return prey
end


# calculate number of predator species
function predspec(bodymass, nP, optimalmassratio, gamma)
	pred = Vector{Float64}()

	for i in 1:length(bodymass)
		tmp = length(findall(feedingeff(bodymass, nP, optimalmassratio, gamma)[i, :] .!= zeros()[]))
		push!(pred, tmp)
	end

	return pred
end


# determine species' diets
function diettype(bodymass, nP, optimalmassratio, gamma)

	linkstructure = feedingeff(bodymass, nP, optimalmassratio, gamma)

	producer = [repeat([false], outer = length(bodymass)-nP)... , repeat([true], outer = nP)...]
	herbivorous = map(x -> sum(linkstructure[producer, x]) > 0, 1:length(bodymass))
	carnivorous = map(x -> sum(linkstructure[.!producer, x]) > 0, 1:length(bodymass))

	omni  = carnivorous .* herbivorous
	carni = carnivorous .* .!omni
	herbi = herbivorous .* .!omni

	diettype = Vector{String}(undef, length(bodymass))

	if sum((producer .+ herbi .+ carni .+ omni) .!= 1) > 0
		println("Overlap between trophic groups or unlinked species!")
	else
		diettype[producer] .= "producer"
		diettype[herbi] .= "herbi"
		diettype[omni]  .= "omni"
		diettype[carni] .= "carni"
	end

	return diettype
end


# calculate capture coefficient based on diet
function capturecoef(bodymass, nP, optimalmassratio, gamma)

	diettype_tmp = diettype(bodymass, nP, optimalmassratio, gamma) # diets: carni, omni, herbi, producer

	beta = Vector{Float64}(undef, length(bodymass))
	beta[diettype_tmp .== "producer"] .= zeros()
	beta[diettype_tmp .== "carni"] = sigmalimit(0.42, 0.05, n = sum(diettype_tmp .== "carni"))
	beta[diettype_tmp .== "omni"]  = sigmalimit(0.19, 0.04, n = sum(diettype_tmp .== "omni"))
	beta[diettype_tmp .== "herbi"] = sigmalimit(0.19, 0.04, n = sum(diettype_tmp .== "herbi"))

	b0 = Vector{Float64}(undef, length(bodymass))
	b0[diettype_tmp .== "producer"] .= zeros()
	b0[diettype_tmp .== "carni"] .= 50
	b0[diettype_tmp .== "omni"]  .= 100
	b0[diettype_tmp .== "herbi"] .= 200

	# calculate capture coefficient
	capturecoef = bodymass.^beta * (b0 .* bodymass.^beta)' .* feedingeff(bodymass, nP, optimalmassratio, gamma)
end


# calculate handling time
function handling(bodymass)
	h_0 = 0.4
	eta_i = sigmalimit(-0.48, 0.03, n = length(bodymass))
	eta_j = sigmalimit(-0.66, 0.02, n = length(bodymass))

	handling = hcat(map(x -> h_0 .* bodymass[x].^eta_i[x] .* bodymass.^eta_j, 1:length(bodymass))...)
end


# calculate feedingrate
function feedingrate(u, omega_i, b_ij, q, int, h_i, bodymass)
	hcat(map(x -> ((omega_i[x] .* b_ij[:, x] .* u.^(1 .+q[:, x])) ./
					(1 .+ int[x] .* (u[x]) .+ omega_i[x] .* h_i[:, x] .* sum(b_ij[:, x] .* u.^(1 .+q[:, x]))) ./ bodymass[x]),
				1:length(bodymass))...)
end
