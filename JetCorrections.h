#ifndef JETCORRECTIONS_H
#define JETCORRECTIONS_H

#include <iostream>
#include <string>
#include <cmath>
#include <TRandom3.h>
#include <TLorentzVector.h>

class JetCorrections {
public:
    JetCorrections(const std::string& jercFile="",
                   const std::string& jetIDFile="",
                   const std::string& vetoMapFile="") {
        std::cout << "JetCorrections initialized." << std::endl;
    }

    double getCorrectedJetPt(double rawPt, double eta, double phi, double rawFactor, bool isMC) {
        // Apply simple raw factor correction
        double correctedPt = rawPt * (1.0 - rawFactor);

        // Add a small MC smearing if needed
        if (isMC) {
            double sigma = 0.02; // 2% smearing (placeholder)
            correctedPt *= gRandom->Gaus(1.0, sigma);
        }

        return correctedPt;
    }

    bool passJetID(const TLorentzVector& jet) {
        return jet.Pt() > 20.0 && std::fabs(jet.Eta()) < 5.0;
    }
};

#endif
