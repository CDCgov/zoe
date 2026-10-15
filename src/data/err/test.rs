use crate::data::err::{DisplayErrStack, ErrorWithContext, WithErrorContext, WithSubitem};

#[test]
fn description_only() {
    let error = ErrorWithContext::new("Failed to parse integer");
    assert_eq!(error.to_string(), "Failed to parse integer");
    assert_eq!(error.display_stack().to_string(), "  → Failed to parse integer\n");
}

#[test]
fn description_and_subitem() {
    let error = ErrorWithContext::new("Failed to parse integer").with_subitem("This is a subitem");
    assert_eq!(error.to_string(), "Failed to parse integer\n| This is a subitem");
    assert_eq!(
        error.display_stack().to_string(),
        "  → Failed to parse integer\n    | This is a subitem\n"
    );
}

#[test]
fn description_and_2subitem_implicit() {
    let error = ErrorWithContext::new("Failed to parse integer").with_subitem("This is a subitem\nThis is another subitem");
    assert_eq!(
        error.to_string(),
        "Failed to parse integer\n| This is a subitem\n| This is another subitem"
    );
    assert_eq!(
        error.display_stack().to_string(),
        "  → Failed to parse integer\n    | This is a subitem\n    | This is another subitem\n"
    );
}

#[test]
fn description_and_2subitem_explicit() {
    let error = ErrorWithContext::new("Failed to parse integer")
        .with_subitem("This is a subitem")
        .with_subitem("This is another subitem");
    assert_eq!(
        error.to_string(),
        "Failed to parse integer\n| This is a subitem\n| This is another subitem"
    );
    assert_eq!(
        error.display_stack().to_string(),
        "  → Failed to parse integer\n    | This is a subitem\n    | This is another subitem\n"
    );
}

#[test]
fn description_and_source() {
    let error = "a".parse::<usize>().unwrap_err().with_context("Failed to parse integer");
    assert_eq!(error.to_string(), "Failed to parse integer");
    assert_eq!(
        error.display_stack().to_string(),
        "  → Failed to parse integer\n  → invalid digit found in string\n"
    );
}

#[test]
fn description_source_and_subitem() {
    let error = "a"
        .parse::<usize>()
        .unwrap_err()
        .with_context("Failed to parse integer")
        .with_subitem("This is a subitem");
    assert_eq!(error.to_string(), "Failed to parse integer\n| This is a subitem");
    assert_eq!(
        error.display_stack().to_string(),
        "  → Failed to parse integer\n    | This is a subitem\n  → invalid digit found in string\n"
    );
}

#[test]
fn description_source_and_2subitem_implicit() {
    let error = "a"
        .parse::<usize>()
        .unwrap_err()
        .with_context("Failed to parse integer")
        .with_subitem("This is a subitem\nThis is another subitem");
    assert_eq!(
        error.to_string(),
        "Failed to parse integer\n| This is a subitem\n| This is another subitem"
    );
    assert_eq!(
        error.display_stack().to_string(),
        "  → Failed to parse integer\n    | This is a subitem\n    | This is another subitem\n  → invalid digit found in string\n"
    );
}

#[test]
fn description_source_and_2subitem_explicit() {
    let error = "a"
        .parse::<usize>()
        .unwrap_err()
        .with_context("Failed to parse integer")
        .with_subitem("This is a subitem")
        .with_subitem("This is another subitem");
    assert_eq!(
        error.to_string(),
        "Failed to parse integer\n| This is a subitem\n| This is another subitem"
    );
    assert_eq!(
        error.display_stack().to_string(),
        "  → Failed to parse integer\n    | This is a subitem\n    | This is another subitem\n  → invalid digit found in string\n"
    );
}

#[test]
fn wrapped_only() {
    let error = ErrorWithContext::from(std::io::Error::other("a".parse::<usize>().unwrap_err()));
    assert_eq!(error.to_string(), "invalid digit found in string");
    assert_eq!(error.display_stack().to_string(), "  → invalid digit found in string\n");
}

#[test]
fn wrapped_and_subitem() {
    let error =
        ErrorWithContext::from(std::io::Error::other("a".parse::<usize>().unwrap_err())).with_subitem("This is a subitem");
    assert_eq!(error.to_string(), "invalid digit found in string\n| This is a subitem");
    assert_eq!(
        error.display_stack().to_string(),
        "  → invalid digit found in string\n    | This is a subitem\n"
    );
}

#[test]
fn wrapped_and_2subitems_implicit() {
    let error = ErrorWithContext::from(std::io::Error::other("a".parse::<usize>().unwrap_err()))
        .with_subitem("This is a subitem\nThis is another subitem");
    assert_eq!(
        error.to_string(),
        "invalid digit found in string\n| This is a subitem\n| This is another subitem"
    );
    assert_eq!(
        error.display_stack().to_string(),
        "  → invalid digit found in string\n    | This is a subitem\n    | This is another subitem\n"
    );
}

#[test]
fn wrapped_and_2subitems_explicit() {
    let error = ErrorWithContext::from(std::io::Error::other("a".parse::<usize>().unwrap_err()))
        .with_subitem("This is a subitem")
        .with_subitem("This is another subitem");
    assert_eq!(
        error.to_string(),
        "invalid digit found in string\n| This is a subitem\n| This is another subitem"
    );
    assert_eq!(
        error.display_stack().to_string(),
        "  → invalid digit found in string\n    | This is a subitem\n    | This is another subitem\n"
    );
}

#[test]
fn wrapped_io() {
    let error = std::io::Error::other("an io error occurred");
    let error = ErrorWithContext::from(error);
    let error = std::io::Error::from(error);
    let error = ErrorWithContext::from(error);
    let error = std::io::Error::from(error);

    assert_eq!(error.to_string(), "an io error occurred");
    assert_eq!(error.display_stack().to_string(), "  → an io error occurred\n");
}

#[test]
fn wrapped_io_code() {
    let error = std::io::Error::from_raw_os_error(22);
    let error = ErrorWithContext::from(error);
    let error = std::io::Error::from(error);
    let error = ErrorWithContext::from(error);
    let error = std::io::Error::from(error);

    assert_eq!(error.raw_os_error(), Some(22));
}

#[test]
fn wrapped_context() {
    let error = ErrorWithContext::new("an io error occurred");
    let error = std::io::Error::from(error);
    let error = ErrorWithContext::from(error);
    let error = std::io::Error::from(error);
    let error = ErrorWithContext::from(error);

    assert_eq!(error.to_string(), "an io error occurred");
    assert_eq!(error.display_stack().to_string(), "  → an io error occurred\n");
}

#[test]
fn more_context() {
    let error = std::io::Error::from(
        "a".parse::<usize>()
            .unwrap_err()
            .with_context("Failed to parse integer")
            .with_subitem("This is a subitem")
            .with_subitem("This is another subitem")
            .with_context("Failed to parse row"),
    )
    .with_context("Failed to parse section")
    .with_subitem("A subitem about the section")
    .with_context("Failed to parse file");

    assert_eq!(error.to_string(), "Failed to parse file");
    assert_eq!(
        error.display_stack().to_string(),
        "  → Failed to parse file\n  → Failed to parse section\n    | A subitem about the section\n  → Failed to parse row\n  → Failed to parse integer\n    | This is a subitem\n    | This is another subitem\n  → invalid digit found in string\n"
    );
}
