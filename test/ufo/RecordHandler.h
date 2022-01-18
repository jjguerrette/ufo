/*
 * (C) Crown copyright 2021, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_UFO_RECORDHANDLER_H_
#define TEST_UFO_RECORDHANDLER_H_

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"

#include "ioda/ObsSpace.h"

#include "oops/runs/Test.h"
#include "oops/util/Expect.h"

#include "ufo/utils/RecordHandler.h"

namespace ufo {
namespace test {

void testRecordHandler(const eckit::LocalConfiguration &conf) {
  util::DateTime bgn(conf.getString("window begin"));
  util::DateTime end(conf.getString("window end"));

  const eckit::LocalConfiguration obsSpaceConf(conf, "obs space");
  ioda::ObsTopLevelParameters obsParams;
  obsParams.validateAndDeserialize(obsSpaceConf);
  ioda::ObsSpace obsspace(obsParams, oops::mpi::world(), bgn, end, oops::mpi::myself());

  // Create a record handler.
  ufo::RecordHandler recordHandler(obsspace);

  // (1) Get input vector.
  // (There is not a getBoolVector option for eckit::Configuration, which is why the conversion from
  // int to bool is performed.)
  const std::vector<int> inputInt = conf.getIntVector("input");
  const std::vector<bool> input(inputInt.begin(), inputInt.end());

  // (2) Check treatment of apply vector

  // Get modified apply vector.
  const std::vector<bool> modified_apply = recordHandler.changeApplyIfRecordsAreSingleObs(input);

  // Get expected modified apply vector.
  const std::vector<int> expected_applyInt = conf.getIntVector("expected_apply");
  const std::vector<bool> expected_apply(expected_applyInt.begin(), expected_applyInt.end());

  EXPECT(expected_apply.size() == modified_apply.size());
  for (size_t idx = 0; idx < modified_apply.size(); ++idx)
    EXPECT(modified_apply[idx] == expected_apply[idx]);

  // (3) Check treatment of isThinned vector.

  // Get modified isThinned vector.
  const std::vector<bool> modified_isThinned =
    recordHandler.changeThinnedIfRecordsAreSingleObs(input);

  // Get expected modified isThinned vector.
  const std::vector<int> expected_isThinnedInt = conf.getIntVector("expected_isThinned");
  const std::vector<bool> expected_isThinned(expected_isThinnedInt.begin(),
                                             expected_isThinnedInt.end());

  EXPECT(expected_isThinned.size() == modified_isThinned.size());
  for (size_t idx = 0; idx < modified_isThinned.size(); ++idx)
    EXPECT(modified_isThinned[idx] == expected_isThinned[idx]);
}

class RecordHandler : public oops::Test {
 private:
  std::string testid() const override {return "ufo::test::RecordHandler";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    const eckit::LocalConfiguration conf(::test::TestEnvironment::config());
    for (const std::string & testCaseName : conf.keys())
    {
      const eckit::LocalConfiguration testCaseConf(::test::TestEnvironment::config(), testCaseName);
      ts.emplace_back(CASE("ufo/RecordHandler/" + testCaseName, testCaseConf)
                      {
                        testRecordHandler(testCaseConf);
                      });
    }
  }

  void clear() const override {}
};

}  // namespace test
}  // namespace ufo

#endif  // TEST_UFO_RECORDHANDLER_H_

