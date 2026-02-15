#include "wheel.h"

#include<algorithm>
#include<bit>
#include<cstring>
#include<numeric>

namespace calcprime{
namespace{

constexpr std::size_t kPresieveBlockWords=16;

std::vector<std::uint16_t> build_presieve_primes(WheelType type){
	(void)type;
	// Keep the presieve modulus compact so all wheels can share fast tables
	// while pre-sieving more than the base wheel factors.
	return {3,5,7,11,13};
}

SmallPrimePattern build_small_pattern(std::uint32_t prime){
	SmallPrimePattern pattern{};
	pattern.prime=prime;
	pattern.phase_count=prime;
	pattern.masks.resize(prime);
	pattern.next_phase.resize(prime);
	std::uint32_t word_stride=static_cast<std::uint32_t>(128%prime);
	pattern.word_stride=word_stride;
	std::uint32_t inv2=(prime+1)/2;
	for(std::size_t bit=0;bit<pattern.start_phase.size();++bit){
		std::uint32_t twice=static_cast<std::uint32_t>(bit<<1);
		while(twice>=prime){
			twice-=prime;
		}
		std::uint32_t phase=(prime-twice);
		if(phase==prime){
			phase=0;
		}
		pattern.start_phase[bit]=static_cast<std::uint8_t>(phase);
	}
	for(std::uint32_t residue=0;residue<prime;++residue){
		std::uint64_t mask=0;
		std::uint32_t offset=static_cast<std::uint32_t>(
			((prime-residue)%prime)*static_cast<std::uint64_t>(inv2)%prime);
		while(offset<64){
			mask|=(1ULL<<offset);
			offset+=prime;
		}
		pattern.masks[residue]=mask;
		pattern.next_phase[residue]=
			static_cast<std::uint32_t>((residue+word_stride)%prime);
	}
	return pattern;
}

Wheel build_wheel(std::uint32_t modulus,WheelType type){
	Wheel wheel;
	wheel.type=type;
	wheel.modulus=modulus;
	wheel.presieved_primes=build_presieve_primes(type);
	wheel.presieve_modulus=1;
	for(std::uint16_t prime : wheel.presieved_primes){
		wheel.presieve_modulus*=prime;
	}
	wheel.allowed.assign(modulus,0);

	for(std::uint32_t r=0;r<modulus;++r){
		if(std::gcd(r,modulus)==1){
			wheel.allowed[r]=1;
			wheel.residues.push_back(static_cast<std::uint16_t>(r));
		}
	}
	if(!wheel.residues.empty()){
		wheel.steps.reserve(wheel.residues.size());
		for(std::size_t i=0;i<wheel.residues.size();++i){
			std::uint32_t a=wheel.residues[i];
			std::uint32_t b=wheel.residues[(i+1)%wheel.residues.size()];
			std::uint32_t step=(b+modulus-a)%modulus;
			if(step==0){
				step=modulus;
			}
			wheel.steps.push_back(static_cast<std::uint16_t>(step));
		}
	}

	static const std::uint32_t kSmallPrimes[]={3, 5, 7, 11,13,17,19,
											   23,29,31,37,41,43,47};
	std::uint32_t small_limit=19u;
	switch(type){
	case WheelType::Mod30:
		small_limit=19u;
		break;
	case WheelType::Mod210:
		small_limit=47u;
		break;
	case WheelType::Mod1155:
		small_limit=47u;
		break;
	}
	for(std::uint32_t prime : kSmallPrimes){
		if(prime>small_limit){
			break;
		}
		wheel.small_patterns.push_back(build_small_pattern(prime));
	}

	std::uint32_t presieve_modulus=wheel.presieve_modulus;
	std::vector<std::uint8_t> presieve_allowed(presieve_modulus,0);
	for(std::uint32_t residue=0;residue<presieve_modulus;++residue){
		if(std::gcd(residue,presieve_modulus)==1){
			presieve_allowed[residue]=1;
		}
	}

	wheel.presieve_word_masks.assign(presieve_modulus,0);
	wheel.presieve_next_phase.assign(presieve_modulus,0);
	for(std::uint32_t phase=0;phase<presieve_modulus;++phase){
		std::uint32_t rem=phase;
		std::uint64_t mask=0;
		for(std::uint32_t bit=0;bit<64;++bit){
			if(!presieve_allowed[rem]){
				mask|=(1ULL<<bit);
			}
			rem+=2;
			if(rem>=presieve_modulus){
				rem-=presieve_modulus;
			}
		}
		wheel.presieve_word_masks[phase]=mask;
		wheel.presieve_next_phase[phase]=static_cast<std::uint16_t>(rem);
	}
	wheel.presieve_block_masks.assign(
		static_cast<std::size_t>(presieve_modulus)*kPresieveBlockWords,0);
	wheel.presieve_next_block_phase.assign(presieve_modulus,0);
	for(std::uint32_t phase=0;phase<presieve_modulus;++phase){
		std::uint32_t rem=phase;
		std::size_t base=
			static_cast<std::size_t>(phase)*kPresieveBlockWords;
		for(std::size_t word=0;word<kPresieveBlockWords;++word){
			wheel.presieve_block_masks[base+word]=wheel.presieve_word_masks[rem];
			rem=wheel.presieve_next_phase[rem];
		}
		wheel.presieve_next_block_phase[phase]=
			static_cast<std::uint16_t>(rem);
	}

	return wheel;
}

} // namespace

const Wheel&get_wheel(WheelType type){
	static const Wheel wheel30=build_wheel(30,WheelType::Mod30);
	static const Wheel wheel210=build_wheel(210,WheelType::Mod210);
	static const Wheel wheel1155=build_wheel(1155,WheelType::Mod1155);
	switch(type){
	case WheelType::Mod30:
		return wheel30;
	case WheelType::Mod210:
		return wheel210;
	case WheelType::Mod1155:
		return wheel1155;
	}
	return wheel30;
}

void Wheel::apply_presieve(std::uint64_t start_value,std::size_t bit_count,
						   std::uint64_t*bits) const{
	if(presieve_modulus==0||presieve_word_masks.empty()||
	   presieve_next_phase.empty()){
		return;
	}
	std::uint32_t phase=
		static_cast<std::uint32_t>(start_value%presieve_modulus);
	std::size_t full_words=bit_count/64;
	std::size_t rem_bits=bit_count%64;
	std::size_t word=0;
	if(!presieve_block_masks.empty()&&!presieve_next_block_phase.empty()){
		std::size_t block_count=full_words/kPresieveBlockWords;
		std::uint64_t*dst=bits;
		for(std::size_t block=0;block<block_count;++block){
			const std::uint64_t*masks=
				&presieve_block_masks[static_cast<std::size_t>(phase)*
									 kPresieveBlockWords];
			for(std::size_t i=0;i<kPresieveBlockWords;++i){
				dst[i]|=masks[i];
			}
			dst+=kPresieveBlockWords;
			phase=presieve_next_block_phase[phase];
		}
		word=block_count*kPresieveBlockWords;
	}
	for(;word<full_words;++word){
		bits[word]|=presieve_word_masks[phase];
		phase=presieve_next_phase[phase];
	}
	if(rem_bits){
		std::uint64_t tail_mask=(1ULL<<rem_bits)-1ULL;
		bits[full_words]|=(presieve_word_masks[phase]&tail_mask);
	}
}

void Wheel::fill_presieve(std::uint64_t start_value,std::size_t bit_count,
						  std::uint64_t*bits) const{
	if(presieve_modulus==0||presieve_word_masks.empty()||
	   presieve_next_phase.empty()){
		return;
	}
	std::uint32_t phase=
		static_cast<std::uint32_t>(start_value%presieve_modulus);
	std::size_t full_words=bit_count/64;
	std::size_t rem_bits=bit_count%64;
	std::size_t word=0;
	if(!presieve_block_masks.empty()&&!presieve_next_block_phase.empty()){
		std::size_t block_count=full_words/kPresieveBlockWords;
		std::uint64_t*dst=bits;
		for(std::size_t block=0;block<block_count;++block){
			const std::uint64_t*masks=
				&presieve_block_masks[static_cast<std::size_t>(phase)*
									 kPresieveBlockWords];
			std::memcpy(dst,masks,sizeof(std::uint64_t)*kPresieveBlockWords);
			dst+=kPresieveBlockWords;
			phase=presieve_next_block_phase[phase];
		}
		word=block_count*kPresieveBlockWords;
	}
	for(;word<full_words;++word){
		bits[word]=presieve_word_masks[phase];
		phase=presieve_next_phase[phase];
	}
	if(rem_bits){
		std::uint64_t tail_mask=(1ULL<<rem_bits)-1ULL;
		bits[full_words]=(presieve_word_masks[phase]&tail_mask);
	}
}

} // namespace calcprime










