#include "State.h"

State::State()
{
	this->isSoft = false;
	this->type = 0;
	this->value = -100.0;
}

State::State (int type) 
{
    this->isSoft = false;
    this->type = type;

    // 0-normal, 1-bust, 2-push, 3-win, 4-lose, 5-blackjack
	if (type == 0){
    	this->value = -100.0;
    }
    if (type == 1){
    	this->value = -1.0;
    }
    if (type == 2){
    	this->value = 0.0;
    }
    if (type == 3){
    	this->value = 1.0;
    }
    if (type == 4){
    	this->value = -1.0;
    }
    if (type == 5){
    	this->value = 1.5;
    }
}

State::State (int player_val, int dealer_upCard)
{
	this->player_val = player_val;
	this->dealer_upCard = dealer_upCard;
	this->type = 0;
	this->isSoft = false;
	this->value = -100.0;
}
