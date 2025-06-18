# For specifying paths

# Juliet's path
# path_to_box <- "~/Library/CloudStorage/Box-Box/"

# # Joey's path
# path_to_box <- "C:/Users/Joseph Fong/Box/"

platform <- sessionInfo()$platform
# Yingyan's path
if (str_detect(platform, "aarch64-apple-darwin20")) {
  path_to_box <- "/Users/yingyanwu/Library/CloudStorage/Box-Box/"
  } else if (str_detect(platform, "x86_64-w64-mingw32/x64")) {
    path_to_box <- "C:/Users/Yingyan Wu/Box/"
  }
