#' @importFrom basilisk BasiliskEnvironment
trex_env <- BasiliskEnvironment(
    envname = "trex_env",
    pkgname = "Trex",
    packages = c(
        "numpy==1.26.4",
        "scipy==1.14.1",
        "h5py==3.12.1"
    ),
    pip = c(
        "tensorflow==2.18.0",
        "keras==3.6.0"
    )
)
