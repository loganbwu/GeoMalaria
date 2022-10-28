test_that("Dynamic data frame class works", {
  df = DynamicFrame$new()
  
  # test that df has the DynamicFrame class
  expect_s3_class(df, "DynamicFrame")
  
  # test that df has no rows of data
  expect_identical(nrow(df$df), 0L)
  
  # initialise
  df2 = DynamicFrame$new(mtcars)
  
  # test that initialising with a dataframe works
  expect_identical(names(df2$df), names(mtcars))
  expect_identical(nrow(df2$df), nrow(mtcars))
  
  # test that initialising with raw data works
  df3 = DynamicFrame$new(x = 0, y = "hi")
  expect_equal(df3$df$x, 0)
  expect_equal(df3$df$y, "hi")
  
  # test appending data
  df3$append(x = 1, y = "bye")
  expect_equal(df3$df$x[2], 1)
  expect_equal(df3$df$y[2], "bye")
})
