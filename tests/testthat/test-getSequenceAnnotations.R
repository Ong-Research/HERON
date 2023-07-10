test_that("getSequenceAnnotations works", {
    probe_meta <- data.frame(
        PROBE_ID=c("A;1","A;2"),
        PROBE_SEQUENCE = c("MSGSASFEGGVFSPYL","SGSASFEGGVFSPYLT")
    )
    expect_equal(
        getSequenceAnnotations("A_1_2", probe_meta),
        data.frame(
            row.names = c("A_1_2"),
            OverlapSeqLength = 15,
            FullSeqStart = 1,
            FullSeqStop = 17,
            FullSeqLength = 17,
            FirstSeq = "MSGSASFEGGVFSPYL",
            LastSeq = "SGSASFEGGVFSPYLT",
            OverlapSeq = "SGSASFEGGVFSPYL",
            FullSeq = "MSGSASFEGGVFSPYLT"
        )
    )
    expect_equal(
        getSequenceAnnotations("A_1_1", probe_meta),
        data.frame(
            row.names = c("A_1_1"),
            OverlapSeqLength = 16,
            FullSeqStart = 1,
            FullSeqStop = 16,
            FullSeqLength = 16,
            FirstSeq = "MSGSASFEGGVFSPYL",
            LastSeq = "MSGSASFEGGVFSPYL",
            OverlapSeq = "MSGSASFEGGVFSPYL",
            FullSeq = "MSGSASFEGGVFSPYL"
        )
    )

    probe_meta2 <- data.frame(
        PROBE_ID=c("A;1", "A;2", "A;3", "A;4"),
        PROBE_SEQUENCE = c("MS","SG", "GS", "SA")
    )

    expect_equal(
        getSequenceAnnotations("A_1_4", probe_meta2),
        data.frame(
            row.names = c("A_1_4"),
            OverlapSeqLength = 0,
            FullSeqStart = 1,
            FullSeqStop = 5,
            FullSeqLength = 5,
            FirstSeq = "MS",
            LastSeq = "SA",
            OverlapSeq = "",
            FullSeq = "MSGSA"

        )
    )

})
