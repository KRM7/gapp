/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GAPP_CROSSOVER_REAL_HPP
#define GAPP_CROSSOVER_REAL_HPP

#include "crossover_base.decl.hpp"
#include "../core/candidate.hpp"
#include "../encoding/gene_types.hpp"
#include "../utility/bounded_value.hpp"

/** Predefined crossover operators for the real encoded genetic algorithm. */
namespace gapp::crossover::real
{
    /**
    * %Arithmetic crossover operator for the real-encoded %GA.
    * 
    * The children are the linear combinations of the parents, such that:
    *   \f[ child_1 = \beta p_1 + (1 - \beta) p_2 \f] 
    *   \f[ child_2 = (1 - \beta) p_1 + \beta p_2 \f]
    * where \f$ \beta \f$ is a random number generated from a uniform distribution on [0.0, 1.0].
    * The same \f$ \beta \f$ value is used for each pair of parent genes.
    */
    class Arithmetic final : public Crossover<RealGene>
    {
    public:
        using Crossover::Crossover;
    private:
        CandidatePair<GeneType> crossover(const GaInfo& ga, const Candidate<GeneType>& parent1, const Candidate<GeneType>& parent2) const override;
    };

    /**
    * BLX-\f$ \alpha \f$ (blend) crossover operator for the real-encoded %GA.
    * 
    * The genes of the children are chosen randomly from a uniform distribution based on
    * the values of the corresponding genes of the parents.
    * 
    * The intervals the child genes are chosen from are:
    *   \f[ [-\alpha I + \min(p_1, p_2),\ \max(p_1, p_2) + \alpha I]
    *       \textrm{,  where } I = |p_1 - p_2| \f]
    * 
    * This crossover operator has a single parameter (alpha), which controls the length
    * of the intervals the child genes are chosen from. Larger alpha values correspond to 
    * larger intervals. The recommended value of alpha is around 0.5.
    */
    class BLXa final : public Crossover<RealGene>
    {
    public:
        /** Create a BLX-alpha crossover operator. */
        constexpr BLXa() noexcept :
            alpha_(GeneType{ 0.5 })
        {}

        /**
        * Create a BLX-alpha crossover operator with the specified parameters.
        *
        * @param pc The crossover probability.
        * @param alpha The alpha parameter of the crossover. Must be a finite non-negative value.
        */
        constexpr explicit BLXa(Probability pc, NonNegative<GeneType> alpha = 0.5) noexcept :
            Crossover(pc), alpha_(alpha)
        {}

        /**
        * Set the alpha parameter of the crossover.
        * This parameter controls the length of the interval the children's genes are
        * chosen from. Larger values values correspond to larger intervals.
        *
        * @param alpha The alpha parameter of the crossover. Must be a finite non-negative value.
        */
        constexpr void alpha(NonNegative<GeneType> alpha) noexcept { alpha_ = alpha; }

        /** @returns The value of the alpha parameter. */
        [[nodiscard]]
        constexpr GeneType alpha() const noexcept { return alpha_; }

    private:
        CandidatePair<GeneType> crossover(const GaInfo& ga, const Candidate<GeneType>& parent1, const Candidate<GeneType>& parent2) const override;

        NonNegative<GeneType> alpha_;
    };

    /**
    * Simulated binary crossover (SBX) operator for the real-encoded %GA.
    * This crossover operator is based on the single-point crossover used in the binary encoded %GA.
    * 
    * This crossover has a single parameter (eta), which controls the shape of the probability
    * distribution the child genes are chosen from. Larger eta values lead to children further away
    * from the parents, while smaller values will result in the children being closer to the parents.
    * Typical values for eta are in the range [1.0, 5.0].
    */
    class SimulatedBinary final : public Crossover<RealGene>
    {
    public:
        /** Create a simulated binary crossover operator. */
        constexpr SimulatedBinary() noexcept :
            eta_(GeneType{ 4.0 })
        {}

        /**
        * Create a simulated binary crossover operator with the specified parameters.
        *
        * @param pc The crossover probability.
        * @param eta The shape parameter of the simulated binary crossover.
        *   Must be finite, non-negative value.
        */
        constexpr explicit SimulatedBinary(Probability pc, NonNegative<GeneType> eta = 4.0) noexcept :
            Crossover(pc), eta_(eta)
        {}

        /**
        * Set the shape parameter of the simulated binary crossover.
        * 
        * The value of this parameter controls the shape of the probability distribution
        * the children's genes are chosen from. Larger value will lead to child genes that
        * are more likely to be further away from the parent genes values.
        * Typical values for this parameter are in the range [1.0, 5.0].
        * 
        * @param eta The eta parameter of the crossover. Must be finite, non-negative value.
        */
        constexpr void eta(NonNegative<GeneType> eta) noexcept { eta_ = eta; }

        /** @returns The eta parameter of the operator. */
        [[nodiscard]]
        constexpr GeneType eta() const noexcept { return eta_; }

    private:
        CandidatePair<GeneType> crossover(const GaInfo& ga, const Candidate<GeneType>& parent1, const Candidate<GeneType>& parent2) const override;

        NonNegative<GeneType> eta_;
    };

    /**
    * %Wright's heuristic crossover operator for the real-encoded %GA.
    * 
    * If p1 is the better parent, then the created children are:
    *   \f[ child_1 = p_1 + w_1 * (p_1 - p_2) \f]
    *   \f[ child_2 = p_1 + w_2 * (p_1 - p_2) \f]
    * where \f$ w_1 \f$ and \f$ w_2 \f$ are random weights generated from a uniform
    * distribution on [0.0, 1.0]. The same weights are used for all of the genes.
    */
    class Wright final : public Crossover<RealGene>
    {
    public:
        using Crossover::Crossover;
    private:
        CandidatePair<GeneType> crossover(const GaInfo& ga, const Candidate<GeneType>& parent1, const Candidate<GeneType>& parent2) const override;
    };


} // namespace gapp::crossover::real

#endif // !GAPP_CROSSOVER_REAL_HPP
