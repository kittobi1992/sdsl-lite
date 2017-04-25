#include "sdsl/int_vector.hpp"
#include "sdsl/rl_enc_vector.hpp"
#include "sdsl/util.hpp"
#include "gtest/gtest.h"
#include <vector>
#include <string>
#include <random>
#include <algorithm>

namespace
{

typedef sdsl::rl_enc_vector<>::size_type size_type;
typedef sdsl::rl_enc_vector<>::value_type value_type;

// The fixture for testing class int_vector.
class RlEncVectorTest : public ::testing::Test
{
    protected:

        RlEncVectorTest() : size(100000),
                            random(size),
                            less_runs(size)  {}

        virtual ~RlEncVectorTest() {}

        virtual void SetUp()
        {
            std::mt19937_64 rng;
            {
                std::uniform_int_distribution<uint64_t> distribution(1, 10000000);
                auto dice = bind(distribution, rng);
                for (size_type i=0; i<size; ++i) {
                    random[i] = dice();
                }
            }
            {
                std::uniform_int_distribution<uint64_t> distribution(2, 5);
                auto dice = bind(distribution, rng);
                less_runs[0] = 1;
                for (size_type i=1; i < size; ++i) {
                    if(i % 1000 == 0) {
                        less_runs[i] = less_runs[i-1] + dice();
                    }
                    else {
                        less_runs[i] = less_runs[i-1] + 1;
                    }
                }
            }
        }

        virtual void TearDown() {}

        size_t size;
        sdsl::int_vector<> random;
        sdsl::int_vector<> less_runs;
};

TEST_F(RlEncVectorTest, CheckValuesOfRandomVector) {
    sdsl::rl_enc_vector<> rl_vec(random);
    ASSERT_EQ(rl_vec.size(),random.size());
    for(size_type i = 0; i < size; ++i) {
        ASSERT_EQ(rl_vec[i],random[i]);
    }
}

TEST_F(RlEncVectorTest, CheckValuesOfVectorWithLessRuns) {
    sdsl::rl_enc_vector<> rl_vec(less_runs);
    ASSERT_EQ(rl_vec.size(),less_runs.size());
    for(size_type i = 0; i < size; ++i) {
        ASSERT_EQ(rl_vec[i],less_runs[i]);
    }
}

}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    if (argc < 2) {
        // LCOV_EXCL_START
        std::cout << "Usage: " << argv[0] << " tmp_dir" << std::endl;
        return 1;
        // LCOV_EXCL_STOP
    }
    return RUN_ALL_TESTS();
}