#ifndef SIGNATURE_LAMBDA_CXX
#define SIGNATURE_LAMBDA_CXX

#include <iostream>
#include "SignatureToolBase.h"
#include "VertexToolBase.h"

namespace signature
{

class LambdaSignature : public SignatureToolBase, public VertexToolBase
{
public:
    explicit LambdaSignature(const fhicl::ParameterSet& pset)
        : _MCPproducer{pset.get<art::InputTag>("MCPproducer", "largeant")}
        , _MCTproducer{pset.get<art::InputTag>("MCTproducer", "generator")}
    {
        configure(pset);
    }

    ~LambdaSignature() override = default;
    
    void configure(fhicl::ParameterSet const& pset) override
    {
        SignatureToolBase::configure(pset);
    }

    TVector3 findVertex(art::Event const& evt) const override;

protected:
    void findSignature(art::Event const& evt, Signature& signature, bool& signature_found) override;

private:
    art::InputTag _MCPproducer;
    art::InputTag _MCTproducer;
};

void LambdaSignature::findSignature(art::Event const& evt, Signature& signature, bool& signature_found)
{
    signature.first = SignatureLambda;
   
    auto const &mcp_h = evt.getValidHandle<std::vector<simb::MCParticle>>(_MCPproducer);
    std::map<int, art::Ptr<simb::MCParticle>> mcp_map;
    for (size_t mcp_i = 0; mcp_i < mcp_h->size(); mcp_i++) {
        const art::Ptr<simb::MCParticle> mcp(mcp_h, mcp_i);
        mcp_map[mcp->TrackId()] = mcp;
    }

    auto addDaughterInteractions = [this, &signature, &mcp_map](const art::Ptr<simb::MCParticle>& particle, auto& self) -> void {
        auto daughters = common::GetDaughters(mcp_map.at(particle->TrackId()), mcp_map);
        for (const auto& daugh : daughters) {
            if (daugh->PdgCode() == particle->PdgCode() && this->assessParticle(*daugh)) {
                this->fillSignature(daugh, signature); 
                self(daugh, self); 
            }
        }
    };

    for (const auto &mcp : *mcp_h) 
    {
        if (abs(mcp.PdgCode()) == 3122 && mcp.Process() == "primary" && mcp.EndProcess() == "Decay" /*&& mcp.NumberDaughters() == 2*/ && !signature_found) 
        {
            auto decay = common::GetDaughters(mcp_map.at(mcp.TrackId()), mcp_map);

            std::vector<art::Ptr<simb::MCParticle>> clean_decay;

            for (const auto& elem : decay) 
            {
              if (elem->Process() != "Decay") 
                continue;

              clean_decay.push_back(elem);
            }

            decay = clean_decay;

            if (decay.size() != 2) continue; 

            std::vector<int> exp_decay = {-211, 2212};
            std::vector<int> fnd_decay;
                
            for (const auto &elem : decay) 
                fnd_decay.push_back(elem->PdgCode());

            std::sort(exp_decay.begin(), exp_decay.end());
            std::sort(fnd_decay.begin(), fnd_decay.end());

            if (fnd_decay == exp_decay) 
            {
                bool all_pass = std::all_of(decay.begin(), decay.end(), [&](const auto& elem) {
                    return this->assessParticle(*elem);
                });

                if (all_pass) 
                {
                    signature_found = true;
                    for (const auto &elem : decay) 
                    {
                        const TParticlePDG* info = TDatabasePDG::Instance()->GetParticle(elem->PdgCode());
                        if (info->Charge() != 0.0) 
                        {
                            this->fillSignature(elem, signature);
                            addDaughterInteractions(elem, addDaughterInteractions);
                        }
                    }

                    break;
                }
            }
        }
    }
}

TVector3 LambdaSignature::findVertex(art::Event const& evt) const
{
    auto const &mcp_h = evt.getValidHandle<std::vector<simb::MCParticle>>(_MCPproducer);
    std::map<int, art::Ptr<simb::MCParticle>> mcp_map;
    for (size_t mcp_i = 0; mcp_i < mcp_h->size(); mcp_i++) 
    {
        const art::Ptr<simb::MCParticle> mcp(mcp_h, mcp_i);
        mcp_map[mcp->TrackId()] = mcp;
    }

    for (const auto &mcp : *mcp_h) {
        if (abs(mcp.PdgCode()) == 3122 && mcp.Process() == "primary" && mcp.EndProcess() == "Decay" && mcp.NumberDaughters() == 2) 
        {
            const TLorentzVector& end_position = mcp.EndPosition();
            return TVector3(end_position.X(), end_position.Y(), end_position.Z());
        }
    }

    return TVector3();
}

DEFINE_ART_CLASS_TOOL(LambdaSignature)
} 

#endif
