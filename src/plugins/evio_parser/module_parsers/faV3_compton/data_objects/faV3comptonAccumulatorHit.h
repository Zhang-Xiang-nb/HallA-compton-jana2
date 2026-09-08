#ifndef FAV3COMPTONACCUMULATORHIT_H
#define FAV3COMPTONACCUMULATORHIT_H

/**
 * @class faV3comptonAccumulatorHit 
 * @brief faV3compton accumulator data
 */
class faV3comptonAccumulatorHit {
public:
    uint32_t acc_samp_overflow;
    uint32_t acc_samp_underflow;
    uint32_t acc_type;
    uint32_t acc_chan;
    uint32_t acc_overflow_timestamp1;
    uint32_t acc_overflow_timestamp2;
    uint32_t acc_underflow_timestamp1;
    uint32_t acc_underflow_timestamp2;
    uint32_t acc_sum_nsample1;
    uint32_t acc_sum_nsample2;
    uint32_t acc_sum1;
    uint32_t acc_sum2;
    uint32_t acc_np_nsboverlapped;
    uint32_t acc_np_nonsa;
    uint32_t acc_np_miss;
    // new combined accumulator data members
    uint64_t acc_sum;
    uint64_t acc_sum_nsample;
    /**
     * @brief Default constructor
     * 
     * Initializes all members to zero.
     */
    faV3comptonAccumulatorHit() : acc_samp_overflow(0), acc_samp_underflow(0), acc_type(0), acc_chan(0),
                 acc_overflow_timestamp1(0), acc_overflow_timestamp2(0), acc_underflow_timestamp1(0), acc_underflow_timestamp2(0), 
                 acc_sum_nsample1(0), acc_sum_nsample2(0), acc_sum1(0), acc_sum2(0), acc_np_nsboverlapped(0), acc_np_nonsa(0), acc_np_miss(0), acc_sum(0), acc_sum_nsample(0) {}

};

#endif // FAV3COMPTONACCUMULATORHIT_H

