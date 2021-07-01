def get_from_events(events, branch_name):
    """Return the branch from the events TTree if it exists, else return None."""

    if branch_name in events.fields:
        return events[branch_name]
    else:
        return None


