#ifndef EVENT_LIST_H
#define EVENT_LIST_H

#include <vector>
#include <string>

// Event list class -- If you change this class make sure to re-run the make file
class Event_List {
	public:
		Event_List()=default;
		Event_List(std::string label_, std::string reweighter_) {label = label_; reweighter = reweighter_;}; // constructor
		std::string reweighter;          // Reweighter type e.g. genie, model, reinteraction
		std::string label;               // name of model used e.g. QEMA, ResDecayEta etc.
		// std::string mode;                // type e.g. unisim or multisim
		std::vector<double> N_reweight;  // A vector of reweighted events the size of this will be equal to the number of universes
	
};

#endif
