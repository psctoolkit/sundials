/* -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2024, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------*/

#include <iostream>
#include <string>
#include <sundials/sundials_core.hpp>
#include <nvector/nvector_serial.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#define TTYPE int
#include "stl/sundials_stlvector.h"

class SUNStlVectorPODTest : public testing::Test {
 protected:
  SUNStlVector_int list;

  virtual void SetUp() override {
    list = SUNStlVector_int_New(2);
  }

  virtual void TearDown() override {
    SUNStlVector_int_Destroy(&list);
  }
};

TEST_F(SUNStlVectorPODTest, NewAndDestroy) {
  EXPECT_NE(list, nullptr);
  EXPECT_EQ(list->size, 0);
  EXPECT_EQ(list->capacity, 2);
  SUNStlVector_int_Destroy(&list);
  EXPECT_EQ(list, nullptr);
}

TEST_F(SUNStlVectorPODTest, IsEmpty) {
  EXPECT_TRUE(SUNStlVector_int_IsEmpty(list));
}

TEST_F(SUNStlVectorPODTest, PushBack) {
  int value1 = 10;
  int value2 = 20;

  SUNStlVector_int_PushBack(list, value1);
  EXPECT_EQ(list->size, 1);
  EXPECT_EQ(*SUNStlVector_int_At(list, 0), value1);

  SUNStlVector_int_PushBack(list, value2);
  EXPECT_EQ(list->size, 2);
  EXPECT_EQ(*SUNStlVector_int_At(list, 1), value2);

  // Test resize
  int value3 = 30;
  SUNStlVector_int_PushBack(list, value3);
  EXPECT_EQ(list->size, 3);
  EXPECT_EQ(list->capacity, 3);
  EXPECT_EQ(*SUNStlVector_int_At(list, 2), value3);
}

TEST_F(SUNStlVectorPODTest, At) {
  int value = 10;
  SUNStlVector_int_PushBack(list, value);

  EXPECT_EQ(*SUNStlVector_int_At(list, 0), value);
  // Test out of bounds
  EXPECT_EQ(SUNStlVector_int_At(list, -1), nullptr);
  EXPECT_EQ(SUNStlVector_int_At(list, 2), nullptr);
}

TEST_F(SUNStlVectorPODTest, Set) {
  int value1 = 10;
  int value2 = 20;
  SUNStlVector_int_PushBack(list, value1);
  SUNStlVector_int_PushBack(list, value2);

  EXPECT_EQ(*SUNStlVector_int_At(list, 0), value1);
  SUNStlVector_int_Set(list, 0, value2);
  EXPECT_EQ(*SUNStlVector_int_At(list, 0), value2);
  // Test out of bounds
  SUNStlVector_int_Set(list, -1, value1);  // No effect
  SUNStlVector_int_Set(list, 2, value1);  // No effect
}

TEST_F(SUNStlVectorPODTest, PopBack) {
  int value1 = 10;
  int value2 = 20;
  SUNStlVector_int_PushBack(list, value1);
  SUNStlVector_int_PushBack(list, value2);

  EXPECT_EQ(list->size, 2);
  SUNStlVector_int_PopBack(list);
  EXPECT_EQ(list->size, 1);
  EXPECT_EQ(*SUNStlVector_int_At(list, 0), value1);
  // Pop from empty list
  SUNStlVector_int_PopBack(list);
  EXPECT_EQ(list->size, 0);
}

#define TTYPE N_Vector
#include "stl/sundials_stlvector.h"

class SUNStlVectorComplexTest : public testing::Test {
 protected:
  SUNStlVector_N_Vector list;
  sundials::Context sunctx;

  virtual void SetUp() override {
    list = SUNStlVector_N_Vector_New(2);
  }

  virtual void TearDown() override {
    SUNStlVector_N_Vector_Destroy(&list);
  }
};

TEST_F(SUNStlVectorComplexTest, NewAndDestroy) {
  EXPECT_NE(list, nullptr);
  EXPECT_EQ(list->size, 0);
  EXPECT_EQ(list->capacity, 2);
  SUNStlVector_N_Vector_Destroy(&list);
  EXPECT_EQ(list, nullptr);
}

TEST_F(SUNStlVectorComplexTest, IsEmpty) {
  EXPECT_TRUE(SUNStlVector_N_Vector_IsEmpty(list));
}

TEST_F(SUNStlVectorComplexTest, PushBack) {
  N_Vector value1 = N_VNew_Serial(1, sunctx);
  N_Vector value2 = N_VNew_Serial(2, sunctx);

  SUNStlVector_N_Vector_PushBack(list, value1);
  EXPECT_EQ(list->size, 1);
  EXPECT_EQ(*SUNStlVector_N_Vector_At(list, 0), value1);
  EXPECT_EQ(N_VGetLength(*SUNStlVector_N_Vector_At(list, 0)), N_VGetLength(value1));

  SUNStlVector_N_Vector_PushBack(list, value2);
  EXPECT_EQ(list->size, 2);
  EXPECT_EQ(*SUNStlVector_N_Vector_At(list, 1), value2);
  EXPECT_EQ(N_VGetLength(*SUNStlVector_N_Vector_At(list, 1)), N_VGetLength(value2));

  // Test resize
  N_Vector value3 = N_VNew_Serial(3, sunctx);
  SUNStlVector_N_Vector_PushBack(list, value3);
  EXPECT_EQ(list->size, 3);
  EXPECT_EQ(list->capacity, 3);
  EXPECT_EQ(*SUNStlVector_N_Vector_At(list, 2), value3);
  EXPECT_EQ(N_VGetLength(*SUNStlVector_N_Vector_At(list, 2)), N_VGetLength(value3));

  N_VDestroy(value1);
  N_VDestroy(value2);
  N_VDestroy(value3);
}

TEST_F(SUNStlVectorComplexTest, At) {
  N_Vector value = N_VNew_Serial(1, sunctx);

  SUNStlVector_N_Vector_PushBack(list, value);
  EXPECT_EQ(*SUNStlVector_N_Vector_At(list, 0), value);

  // Test out of bounds
  EXPECT_EQ(SUNStlVector_N_Vector_At(list, -1), nullptr);
  EXPECT_EQ(SUNStlVector_N_Vector_At(list, 2), nullptr);

  N_VDestroy(value);
}

TEST_F(SUNStlVectorComplexTest, Set) {
  N_Vector value1 = N_VNew_Serial(1, sunctx);
  N_Vector value2 = N_VNew_Serial(2, sunctx);

  SUNStlVector_N_Vector_PushBack(list, value1);
  SUNStlVector_N_Vector_PushBack(list, value2);
  EXPECT_EQ(*SUNStlVector_N_Vector_At(list, 0), value1);
  SUNStlVector_N_Vector_Set(list, 0, value2);
  EXPECT_EQ(*SUNStlVector_N_Vector_At(list, 0), value2);

  // Test out of bounds
  SUNStlVector_N_Vector_Set(list, -1, value1);
  SUNStlVector_N_Vector_Set(list, 2, value1);

  N_VDestroy(value1);
  N_VDestroy(value2);
}

TEST_F(SUNStlVectorComplexTest, PopBack) {
  N_Vector value1 = N_VNew_Serial(1, sunctx);
  N_Vector value2 = N_VNew_Serial(2, sunctx);

  SUNStlVector_N_Vector_PushBack(list, value1);
  SUNStlVector_N_Vector_PushBack(list, value2);
  EXPECT_EQ(list->size, 2);

  SUNStlVector_N_Vector_PopBack(list);
  EXPECT_EQ(list->size, 1);
  EXPECT_EQ(*SUNStlVector_N_Vector_At(list, 0), value1);

  // Pop from empty list
  SUNStlVector_N_Vector_PopBack(list);
  EXPECT_EQ(list->size, 0);

  N_VDestroy(value1);
  N_VDestroy(value2);
}
