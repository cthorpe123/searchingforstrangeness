#ifndef SIGNATURE_NEUTRALKAON_CXX
#define SIGNATURE_NEUTRALKAON_CXX

#include <iostream>
#include <optional>
#include "SignatureToolBase.h"
#include "DecayVertexProvider.h"

namespace signature
{

class KaonShortSignature : public SignatureToolBase, public DecayVertexProvider
{
public:
    explicit KaonShortSignature(const fhicl::ParameterSet& pset)
        : _HitProducer{pset.get<art::InputTag>("HitProducer", "gaushit")}
        , _MCPproducer{pset.get<art::InputTag>("MCPproducer", "largeant")}
        , _MCTproducer{pset.get<art::InputTag>("MCTproducer", "generator")}
        , _BacktrackTag{pset.get<art::InputTag>("BacktrackTag", "gaushitTruthMatch")}
    {
        configure(pset);
    }

    ~KaonShortSignature() override = default;

    void configure(fhicl::ParameterSet const& pset) override
    {
        SignatureToolBase::configure(pset);
    }

    TVector3 getDecayVertex(art::Event const& evt) const override;

protected:
    void findSignature(art::Event const& evt, Signature& sig, bool& found_signature) override;

private:
    art::InputTag _HitProducer;
    art::InputTag _MCPproducer;
    art::InputTag _MCTproducer;
    art::InputTag _BacktrackTag;
};

void KaonShortSignature::findSignature(art::Event const& evt, Signature& sig, bool& found_signature)
{
    auto const &mcp_h = evt.getValidHandle<std::vector<simb::MCParticle>>(_MCPproducer);

    std::map<int, art::Ptr<simb::MCParticle>> mcp_map;
    for (size_t d = 0; d < mcp_h->size(); d++) {
        const art::Ptr<simb::MCParticle> mcp(mcp_h, d);
        mcp_map[mcp->TrackId()] = mcp;
    }

    for (const auto &t_part : *mcp_h) {
        if (abs(t_part.PdgCode()) == 311 && t_part.Process() == "primary" && t_part.EndProcess() == "Decay" && t_part.NumberDaughters() == 1 && !found_signature) {
            auto g_dtrs = common::GetDaughters(mcp_map.at(t_part.TrackId()), mcp_map);
            if (g_dtrs.size() != 1) continue; 

            auto g_part = g_dtrs.at(0);
            if (g_part->PdgCode() == 310 && g_part->Process() == "Decay" && g_part->EndProcess() == "Decay" && g_part->NumberDaughters() == 2 && !found_signature)
            {
                auto daughters = common::GetDaughters(mcp_map.at(g_part->TrackId()), mcp_map);
                if (daughters.size() != 2) continue;
                
                std::vector<int> exp_dtrs = {-211, 211};
                std::vector<int> fnd_dtrs;
                
                for (const auto &dtr : daughters) 
                {
                    fnd_dtrs.push_back(dtr->PdgCode());
                    std::cout << dtr->PdgCode() << std::endl;
                }

                std::sort(exp_dtrs.begin(), exp_dtrs.end());
                std::sort(fnd_dtrs.begin(), fnd_dtrs.end());

                if (fnd_dtrs == exp_dtrs) 
                {
                    bool all_above_threshold = std::all_of(daughters.begin(), daughters.end(), [&](const auto& dtr) {
                        return this->aboveThreshold(*dtr);
                    });

                    if (all_above_threshold) 
                    {
                        found_signature = true;
                        for (const auto &dtr : daughters) 
                            this->fillSignature(dtr, sig);

                        break;
                    }
                }
            }
        }
    }
}

TVector3 KaonShortSignature::getDecayVertex(art::Event const& evt) const
{
    auto const &mcp_h = evt.getValidHandle<std::vector<simb::MCParticle>>(_MCPproducer);

    std::map<int, art::Ptr<simb::MCParticle>> mcp_map;
    for (size_t d = 0; d < mcp_h->size(); d++) {
        const art::Ptr<simb::MCParticle> mcp(mcp_h, d);
        mcp_map[mcp->TrackId()] = mcp;
    }

    for (const auto &t_part : *mcp_h) {
        if (abs(t_part.PdgCode()) == 311 && t_part.Process() == "primary" && t_part.EndProcess() == "Decay" && t_part.NumberDaughters() == 1) {
            auto g_dtrs = common::GetDaughters(mcp_map.at(t_part.TrackId()), mcp_map);
            if (g_dtrs.size() != 1) continue; 

            auto g_part = g_dtrs.at(0);
            if (g_part->PdgCode() == 310 && g_part->Process() == "Decay" && g_part->EndProcess() == "Decay" && g_part->NumberDaughters() == 2)
            {
                const TLorentzVector& end_position = g_part->EndPosition();
                return TVector3(end_position.X(), end_position.Y(), end_position.Z());
            }
        }
    }

    return TVector3();
}

DEFINE_ART_CLASS_TOOL(KaonShortSignature)
} 

#endif
